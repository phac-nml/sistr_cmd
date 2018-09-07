import logging
from typing import Dict, Type, IO
import re
import os
import zipfile
from datetime import datetime
from glob import glob

import attr
import pandas as pd


def find_file_in_zip(zip_file: zipfile.ZipFile, pattern: str, escape: bool = False) -> zipfile.ZipInfo:
    regex = re.compile(re.escape(pattern) if escape else pattern)
    for zinfo in zip_file.infolist():
        if not zinfo.is_dir() and regex.search(zinfo.filename):
            return zinfo


def find_file_in_dir(dir_path: str, pattern: str, escape: bool = False) -> str:
    regex = re.compile(re.escape(pattern) if escape else pattern)
    for filename in glob(os.path.join(dir_path, '**'), recursive=True):
        if os.path.isfile(filename) and regex.search(filename):
            return os.path.abspath(filename)


def read_hdf_from_buffer(buffer, key='cgmlst'):
    from pandas import HDFStore
    with HDFStore("data.h5",
                  mode="r",
                  driver="H5FD_CORE",
                  driver_core_backing_store=0,
                  driver_core_image=buffer.read()) as store:
        return store[key]


@attr.s
class SISTRDatabaseMetadata:
    version: str = attr.ib(default=None)
    date_created: datetime = attr.ib(default=None)
    doi: str = attr.ib(default=None)
    url: str = attr.ib(default=None)
    author: str = attr.ib(default=None)
    author_email: str = attr.ib(default=None)


@attr.s
class SISTRDataBase:
    reference_cgmlst_profiles: pd.DataFrame = attr.ib(default=None)
    reference_genome_info: pd.DataFrame = attr.ib(default=None)
    is_zip: bool = attr.ib(default=False)
    db_path: str = attr.ib(default=None)
    metadata: SISTRDatabaseMetadata = attr.ib(default=None)
    _zipfile: zipfile.ZipFile = attr.ib(default=None, repr=False)

    METADATA_FILENAME = 'sistr.yml'
    CGMLST_PROFILES_FILENAME = 'cgmlst-profiles.hdf'
    CGMLST_CENTROID_ALLELES_FILENAME = 'cgmlst-centroid.fasta'
    CGMLST_FULL_ALLELES_FILENAME = 'cgmlst-full.fasta'
    WZX_FILENAME = 'wzx.fasta'
    WZY_FILENAME = 'wzy.fasta'
    FLIC_FILENAME = 'fliC.fasta'
    FLJB_FILENAME = 'fljB.fasta'
    REF_GENOME_FILENAME = 'reference-genomes'
    REF_MASH_SKETCHES_FILENAME = 'sistr.msh'
    SEROTYPE_INFO_TABLE_FILENAME = 'serotype.csv'
    EXPECTED_FILES = [
        CGMLST_PROFILES_FILENAME,
        CGMLST_FULL_ALLELES_FILENAME,
        CGMLST_CENTROID_ALLELES_FILENAME,
        WZX_FILENAME,
        WZY_FILENAME,
        FLIC_FILENAME,
        FLJB_FILENAME,
        REF_GENOME_FILENAME,
        REF_MASH_SKETCHES_FILENAME,
        SEROTYPE_INFO_TABLE_FILENAME,
    ]

    def has_all_expected_files(self):
        if self.is_zip:
            with zipfile.ZipFile(self.db_path) as zf:
                missing_files = [filename for filename in SISTRDataBase.EXPECTED_FILES if
                                 not find_file_in_zip(zf, filename)]
        else:
            missing_files = [filename for filename in SISTRDataBase.EXPECTED_FILES if
                             not find_file_in_dir(self.db_path, filename)]
        if len(missing_files) > 0:
            logging.error(f'SISTR {"zip" if self.is_zip else "directory"} Database "{self.db_path}" is missing '
                          f'"{missing_files}". Cannot run SISTR with this database.')
            raise Exception()

    def _try_parse_metadata(self):
        if self.is_zip:
            with zipfile.ZipFile(self.db_path) as zf:
                metadata_file_info = find_file_in_zip(zf, SISTRDataBase.METADATA_FILENAME)
                if metadata_file_info:
                    with zf.open(metadata_file_info) as zfh:
                        import yaml
                        self.metadata = SISTRDatabaseMetadata(**yaml.load(zfh))
        else:
            metadata_filepath = find_file_in_dir(self.db_path, SISTRDataBase.METADATA_FILENAME)
            if metadata_filepath:
                with open(metadata_filepath) as f:
                    import yaml
                    self.metadata = SISTRDatabaseMetadata(**yaml.load(f))

    def get_file_handle(self, filename_pattern: str) -> IO:
        if self.is_zip:
            if self._zipfile:
                zf = self._zipfile
            else:
                zf = zipfile.ZipFile(self.db_path)
                self._zipfile = zf

            return zf.open(find_file_in_zip(zf, filename_pattern))
        else:
            return open(find_file_in_dir(self.db_path, filename_pattern))

    @classmethod
    def from_zip(cls: Type['SISTRDataBase'], sistr_db_zip_path: str) -> 'SISTRDataBase':
        obj = cls(is_zip=True,
                  db_path=sistr_db_zip_path)
        obj.has_all_expected_files()
        obj._try_parse_metadata()
        obj._zipfile = zipfile.ZipFile(sistr_db_zip_path)
        cgmlst_profiles_zipinfo = find_file_in_zip(obj._zipfile, cls.CGMLST_PROFILES_FILENAME)
        with obj._zipfile.open(cgmlst_profiles_zipinfo) as zfh:
            obj.reference_cgmlst_profiles = read_hdf_from_buffer(zfh)
        ref_genomes_info_fileinfo = find_file_in_zip(obj._zipfile, cls.SEROTYPE_INFO_TABLE_FILENAME)
        with obj._zipfile.open(ref_genomes_info_fileinfo) as zfh:
            if ref_genomes_info_fileinfo.filename.endswith('.csv'):
                obj.reference_genome_info = pd.read_csv(zfh)
            else:
                # try parsing as tab-delimited
                obj.reference_genome_info = pd.read_table(zfh)
        return obj

    @classmethod
    def from_dir(cls: Type['SISTRDataBase'], sistr_db_dir_path: str) -> 'SISTRDataBase':
        obj = cls(is_zip=False,
                  db_path=sistr_db_dir_path)
        obj.has_all_expected_files()
        obj._try_parse_metadata()
        cgmlst_profiles_filepath = find_file_in_dir(obj.db_path, cls.CGMLST_PROFILES_FILENAME)
        obj.reference_cgmlst_profiles = pd.read_hdf(cgmlst_profiles_filepath)
        ref_genomes_info_filepath = find_file_in_dir(obj.db_path, cls.SEROTYPE_INFO_TABLE_FILENAME)
        if ref_genomes_info_filepath.endswith('.csv'):
            obj.reference_genome_info = pd.read_csv(ref_genomes_info_filepath)
        else:
            # try parsing as tab-delimited
            obj.reference_genome_info = pd.read_table(ref_genomes_info_filepath)
        return obj

    def __del__(self):
        if self._zipfile:
            logging.debug(f'Closing ZipFile handle {self._zipfile}')
            self._zipfile.close()
            logging.debug(f'ZipFile handle closed: {self._zipfile}')
