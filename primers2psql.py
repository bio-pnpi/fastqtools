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
import psycopg

from FastQTools.ExpandPrimer import ExpandPrimer

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

def main():
    parser = argparse.ArgumentParser(
        prog="fastq2sql",
        description="load data from fastq to sql"
    )
    parser.add_argument(
        '--primers',
        dest='primersFnm',
        type=str,
        required=True,
        help='File with primers'
    )
    args = parser.parse_args()
    db = psycopg.connect('dbname=nanoseq user=nanoseq')
    cc = db.cursor()

    Primers = readPrimers(args.primersFnm)
    # now add primers to sql
    for primer in Primers:
        try:
            cc.execute("""
                    INSERT INTO dprimers(primer, len)
                    VALUES (%(primer)s, %(len)s)
                    """,
                    {'primer': primer, 'len': len(primer)}
                   )
        except psycopg.errors.UniqueViolation:
            print(f"dPrimer already exist {primer}")
        db.commit()
    print("Primers is: {}".format(Primers))
    for primer in Primers:
        cc.execute("SELECT id FROM dprimers WHERE primer = %(primer)s", {"primer": primer})
        did = cc.fetchone()[0]
        ePrimer = ExpandPrimer(primer)
        ePrimer.Expand()
        # fill table with primers
        for fp in ePrimer.GetForward():
            try:
                cc.execute("""
                        INSERT INTO fprimers(primer, dprimer)
                        VALUES(%(primer)s, %(dprimer)s)""",
                       {'primer': fp, 'dprimer': did})
            except psycopg.errors.UniqueViolation:
                print(f"fPrimer {fp} already exist for dPrimer {primer}")
            db.commit()
        for rp in ePrimer.GetReverse():
            try:
                cc.execute("""
                        INSERT INTO rprimers(primer, dprimer)
                        VALUES(%(primer)s, %(dprimer)s)""",
                       {'primer': rp, 'dprimer': did})
            except psycopg.errors.UniqueViolation:
                print(f"fPrimer {rp} already exist for dPrimer {primer}")
            db.commit()

if __name__ == '__main__':
    main()
