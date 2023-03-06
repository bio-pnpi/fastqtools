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
import psycopg


def checkFastqGz(fnm):
    if fnm.endswith(".fastq.gz"):
        return True
    else:
        return False


def main():
    parser = argparse.ArgumentParser(
        prog="fastq2sql",
        description="load data from fastq to sql"
    )
    parser.add_argument(
        '--fastq',
        type=str,
        nargs='+',
        required=True,
        help='Input file(s) in fastq format'
    )
    args = parser.parse_args()

    db = psycopg.connect('dbname=nanoseq user=nanoseq')
    cc = db.cursor()

    for file in args.fastq:
        fd = gzip.open(file, 'rb')
        data = fd.readlines()
        fd.close()
        # get number of blocks
        nblocks = len(data) // 4
        for block in range(nblocks):
            # parse header
            hdrFull = str(data[4 * block + 0], 'UTF-8').strip().split()
            for h in hdrFull:
                if h.startswith('@'):
                    sid = h[1:]
                elif h.startswith('runid='):
                    runid = h.replace('runid=', '')
                elif h.startswith('sampleid='):
                    sampleid = h.replace('sampleid=', '')
                elif h.startswith('read='):
                    read = int(h.replace('read=', ''))
                elif h.startswith('ch='):
                    ch = int(h.replace('ch=', ''))
                elif h.startswith('start_time='):
                    start_time = h.replace('start_time=', '')
                elif h.startswith('model_version_id='):
                    model_version_id = h.replace('model_version_id=', '')
                elif h.startswith('barcode='):
                    barcode = h.replace('barcode=', '')
                else:
                    print("Error, unknown line")

            cc.execute("""
                    INSERT INTO header(id, runid, sampleid, read, ch, start_time, model_version_id, barcode)
                    VALUES (%(sid)s, %(runid)s, %(sampleid)s, %(read)s, %(ch)s, %(start_time)s, %(model_version_id)s, %(barcode)s)
                    """,
                    {'sid': sid, 'runid': runid, 'sampleid' : sampleid,
                     'read': read, 'ch': ch, 'start_time': start_time,
                     'model_version_id': model_version_id, 'barcode': barcode}
            )
            seq = str(data[4 * block + 1], 'UTF-8').strip()
            quality = str(data[4 * block + 3], 'UTF-8').strip()
            cc.execute("""
                INSERT INTO sequence(id, sequence, quality, len)
                VALUES(%(sid)s, %(seq)s, %(quality)s, %(len)s)
                """,
                {'sid' : sid, 'seq': seq, 'quality': quality, 'len': len(seq) }
            )
            db.commit()
            print('hdr = {}\nseq = {}\nquality = {}'.format(hdrFull, seq, quality))


if __name__ == '__main__':
    main()
