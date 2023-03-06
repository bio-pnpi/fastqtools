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
import gzip

from FastQTools.ExpandPrimer import ExpandPrimer
from FastQTools.StripPrimer import StripPrimer

parser = argparse.ArgumentParser(
    prog="fastq-strip-primers",
    description="Strip primers from fastq.gz files"
)
parser.add_argument(
    '--fastq',
    type=str,
    nargs='+',
    required=True,
    help='Input file(s) in fastq format'
)

parser.add_argument(
    '--primers',
    dest='primersFnm',
    type=str,
    required=True,
    help='File with primers'
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


def main():
    # get primers
    Primers = readPrimers(args.primersFnm)
    print("Primers is: {}".format(Primers))
    fPrimers = []
    rPrimers = []
    for primer in Primers:
        ePrimer = ExpandPrimer(primer)
        ePrimer.Expand()
        for fp in ePrimer.GetForward():
            fPrimers.append(fp)
        for rp in ePrimer.GetReverse():
            rPrimers.append(rp)
    print("Expanded forward primers: {}".format(fPrimers))
    print("Expanded reverse primers: {}".format(rPrimers))
    # process fastq
    for fastq in args.fastq:
        if checkFastqGz(fastq):
            print("Processing {}".format(fastq))
            fd = gzip.open(fastq, 'rb')
            lines = fd.readlines()
            fd.close()
            fd = gzip.open(fastq.replace(".fastq.gz", "-stripped-primers.fastq.gz"), 'wb')
            nblocks = len(lines) // 4
            for block in range(nblocks):
                # read header
                fd.write("{}\n".format(str(lines[4*block+0], 'UTF-8').strip()).encode())
                # read seq
                seq = str(lines[4*block+1], 'UTF-8').strip()
                # read additional
                add = str(lines[4*block+2], 'UTF-8').strip()
                # read quality
                quality = str(lines[4*block+3], 'UTF-8').strip()
                sp = StripPrimer(seq=seq, quality=quality, fprimers=fPrimers, rprimers=rPrimers)
                sp.Strip()
                seq, quality = sp.GetStripped()
                # write back seq
                fd.write("{}\n".format(seq).encode())
                # write back add
                fd.write("{}\n".format(add).encode())
                # write back quality
                fd.write("{}\n".format(quality).encode())
            fd.close()
        else:
            print("Error: input files should be .fastq.gz")
            exit(1)


if __name__ == '__main__':
    main()
