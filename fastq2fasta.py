#!/usr/bin/env python3
#
# This file is part of the FastQTool distribution (https://gitlab.com/bio-pnpi/fastqtools or http://xxx.github.io).
# Copyright (c) 2023 Alexey Shvetsov.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# SPDX-FileCopyrightText: 2023 Alexey Shvetsov <alexxyum@gmail.com>
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse

from FastQTools.ExpandPrimer import ExpandPrimer
from FastQTools.ReadFastqGzToFasta import ReadFasqGzToFasta
from FastQTools.StripPrimer import StripPrimer

parser = argparse.ArgumentParser(
    prog="fastq2fasta",
    description="Read and convert fastq to fasta"
)

parser.add_argument(
    '--fastq',
    type=str,
    nargs='+',
    required=True,
    help='Input file(s) in fastq format'
)

parser.add_argument(
    '--merge',
    type=str,
    help='Merged file with all seq in fasta format'
)

parser.add_argument(
    '--primers',
    dest='primersFnm',
    type=str,
    help='File with primers'
)
parser.add_argument(
    '--strip-primers',
    dest='stripPrimers',
    action='store_true'
)

args = parser.parse_args()


def readPrimers(primersFnm):
    fd = open(primersFnm)
    data = fd.readlines()
    fd.close()
    Primers = []
    for line in data:
        primer = line.strip()
        if checkPrimer(primer):
            Primers.append(line.strip())
        else:
            print("Error: Input primers has malformed lines {}".format(primer))
            exit(1)
    return Primers


def checkPrimer(Primer):
    pAlphabet = frozenset('ATGCMRWSYKVHDBN')
    return pAlphabet.issuperset(Primer)


def checkFastqGz(fnm):
    if fnm.endswith(".fastq.gz"):
        return True
    else:
        return False


def writeFasta(fnm, hdr, seq):
    fo = open(fnm, 'w+')
    num = len(hdr)
    for i in range(num):
        print(">{h}\n{s}".format(h=hdr[i], s=seq[i]), file=fo)
    fo.close()


def main():
    Primers = None
    fPrimers = []
    rPrimers = []
    hdrall = []
    seqall = []
    # deal with primers
    if args.stripPrimers:
        if args.primersFnm is not None:
            Primers = readPrimers(args.primersFnm)
            print("Primers is: {}".format(Primers))
        else:
            print("Error: Strip Primers option is set, while no Primers Given")
            exit(1)
        if Primers is not None:
            # Expand all primers
            for primer in Primers:
                ePrimer = ExpandPrimer(primer)
                ePrimer.Expand()
                for fp in ePrimer.GetForward():
                    fPrimers.append(fp)
                for rp in ePrimer.GetReverse():
                    rPrimers.append(rp)
            print("Expanded forward primers: {}".format(fPrimers))
            print("Expanded reverse primers: {}".format(rPrimers))

    # Now read and convers fastq files
    for fastq in args.fastq:
        if checkFastqGz(fastq):
            print("Processing {}".format(fastq))
            fqgz = ReadFasqGzToFasta(fastq)
            fqgz.Process()
            hdr = fqgz.GetHdrs()
            seq = fqgz.GetSeqs()
            print("Read {} records from {} writing result".format(fqgz.GetNumRecords(), fastq))
            writeFasta(fastq.replace(".fastq.gz", ".fasta"), hdr, seq)
            if args.stripPrimers:
                sp = StripPrimer(seq=seq, quality=None, fprimers=fPrimers, rprimers=rPrimers)
                sp.Strip()
                seq = sp.GetStripped()
            if args.merge is not None:
                for i in range(fqgz.GetNumRecords()):
                    hdrall.append(hdr[i])
                    seqall.append(seq[i])
        else:
            print("Error: input files should be .fastq.gz")
            exit(1)

    if args.merge is not None:
        writeFasta(args.merge, hdrall, seqall)


if __name__ == '__main__':
    main()
