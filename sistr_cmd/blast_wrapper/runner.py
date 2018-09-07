import logging
import os
import re
import shutil
from datetime import datetime
from io import UnsupportedOperation
from subprocess import Popen, PIPE
from typing import IO

import attr
import errno

from sistr_cmd.blast_wrapper.const import BLAST_TABLE_COLS


@attr.s
class BlastRunner:
    fasta_path: str = attr.ib()
    tmp_work_dir: str = attr.ib()
    blast_db_created: bool = attr.ib(default=False)

    def _create_tmp_folder(self) -> None:
        count = 1
        tmp_dir = self.tmp_work_dir
        while True:
            try:
                logging.info('Trying to create analysis directory at: %s', tmp_dir)
                os.makedirs(tmp_dir)
                break
            except OSError as e:
                logging.warning('Error on creation of tmp analysis directory "{}"! {}'.format(
                    tmp_dir,
                    e
                ))

                tmp_dir = '{}_{}'.format(self.tmp_work_dir, count)
                count += 1

        self.tmp_work_dir = tmp_dir

    def _copy_fasta_to_work_dir(self) -> None:
        filename = os.path.basename(self.fasta_path)
        filename_no_spaces = re.sub(r'\W', '_', filename)
        dest_path = os.path.join(self.tmp_work_dir, filename_no_spaces)
        if self.fasta_path == dest_path:
            self.tmp_fasta_path = dest_path
            return
        os.symlink(self.fasta_path, dest_path)
        self.tmp_fasta_path = dest_path

    def _run_makeblastdb(self) -> str:
        work_dir = os.path.dirname(self.tmp_fasta_path)
        filename = os.path.basename(self.tmp_fasta_path)
        nin_filepath = os.path.join(work_dir, filename + '.nin')
        if os.path.exists(nin_filepath):
            self.blast_db_created = True
            return self.tmp_fasta_path

        p = Popen(['makeblastdb',
                   '-in', f'{self.tmp_fasta_path}',
                   '-dbtype', 'nucl'],
                  stdout=PIPE,
                  stderr=PIPE)
        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()
        if stdout is not None and stdout != '':
            logging.debug('makeblastdb on {0} STDOUT: {1}'.format(self.tmp_fasta_path, stdout))
        if stderr is not None and stderr != '':
            logging.debug('makeblastdb on {0} STDERR: {1}'.format(self.tmp_fasta_path, stderr))

        if os.path.exists(nin_filepath):
            self.blast_db_created = True
        else:
            ex_msg = 'makeblastdb was not able to create a BLAST DB for {0}. STDERR: {1}'.format(filename, stderr)
            logging.error(ex_msg)
            raise Exception(ex_msg)

    def filename(self) -> str:
        return os.path.basename(self.fasta_path)

    @staticmethod
    def can_pass_to_stdin(fh):
        try:
            fh.fileno()
            return True
        except UnsupportedOperation:
            return False

    def blast(self, query: IO, blast_task: str = 'megablast', evalue: float = 1e-20,
              min_percent_identity: float = 85.0) -> str:

        if not self.blast_db_created:
            self.prep_blast()

        can_pass_query_to_stdin = BlastRunner.can_pass_to_stdin(query)

        query_filename = os.path.basename(query.name)
        subject_filename = os.path.basename(self.tmp_fasta_path)
        timestamp = f'{datetime.now():%Y%b%d_%H_%M_%S}'
        outfile = os.path.join(self.tmp_work_dir, f'{query_filename}-{subject_filename}-{timestamp}.blast')
        p = Popen(['blastn',
                   '-task', blast_task,
                   '-query', '-',
                   '-db', f'{self.tmp_fasta_path}',
                   '-evalue', f'{evalue}',
                   '-dust', 'no',
                   '-perc_identity', f'{min_percent_identity}',
                   '-out', outfile,
                   '-outfmt', f'6 {" ".join(BLAST_TABLE_COLS)}'],
                  stdout=PIPE,
                  stderr=PIPE,
                  stdin=query if can_pass_query_to_stdin else PIPE)

        if not can_pass_query_to_stdin:
            for l in query:
                try:
                    p.stdin.write(l)
                except IOError as e:
                    if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
                        break
                    else:
                        raise
            p.stdin.close()
        p.wait()

        stdout = p.stdout.read()
        stderr = p.stderr.read()
        if stdout is not None and stdout != '':
            logging.debug(f'blastn on db "{subject_filename}" and query "{query_filename}" STDOUT: {stdout}')
        if stderr is not None and stderr != '':
            logging.debug(f'blastn on db "{subject_filename}" and query "{query_filename}" STDERR: {stderr}')

        if os.path.exists(outfile):
            return outfile
        else:
            ex_msg = (f'blastn output file "{outfile}" was not created! '
                      f'Subject="{subject_filename}", '
                      f'Query="{query_filename}"')
            logging.error(ex_msg)
            raise Exception(ex_msg)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    def cleanup(self):
        if os.path.exists(self.tmp_work_dir):
            shutil.rmtree(self.tmp_work_dir)
        self.blast_db_created = False

    def __enter__(self):
        self.prep_blast()

    def prep_blast(self) -> None:
        self._create_tmp_folder()
        self._copy_fasta_to_work_dir()
        self._run_makeblastdb()

    def run(self, query: IO, **kwargs) -> str:
        self.prep_blast()
        blast_outfile = self.blast(query, **kwargs)
        return blast_outfile
