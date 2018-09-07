#!/bin/sh

set -e

HERE="$PWD"
WGET="wget --quiet"
UNTAR="tar xf"

BLASTVER=2.7.1
BLAST=ncbi-blast-$BLASTVER+
echo "* $BLAST"
$WGET ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLASTVER/$BLAST-x64-linux.tar.gz
$UNTAR $BLAST.tar.bz2
PATH=$HERE/$BLAST/bin:$PATH


MAFFTVER=7.407
MAFFTTGZ=mafft-$MAFFTVER-linux.tgz
echo "* $MAFFTTGZ"
$WGET https://mafft.cbrc.jp/alignment/software/$MAFFTTGZ
$UNTAR $MAFFTTGZ
export MAFFT_BINARIES=$HERE/mafft-linux64/mafftdir/libexec
PATH=$HERE/mafft-linux64/mafftdir/bin:$PATH

MASHVER=2.0
MASH="mash-Linux64-v$MASHVER"
MASHTAR="$MASH.tar"
echo "* $MASH"
$WGET https://github.com/marbl/Mash/releases/download/v$MASHVER/$MASHTAR
$UNTAR $MASHTAR
PATH=$HERE/$MASH:$PATH

# cover anything else
PATH=$PATH:$HERE

echo "Deleting source files"
rm -vf "$HERE/*.tar.*"
rm -vf "$HERE/*.tar"
rm -vf "$HERE/*.tgz"

echo $PATH
export PATH
