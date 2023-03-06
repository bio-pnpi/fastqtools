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

import gzip


class ReadFasqGzToFasta(object):
    def __init__(self, fnm):
        self.fnm = fnm
        self.hdr = []
        self.seq = []
        self.nrecords = 0
        self.lines = None
        self.stag = False

    def readLines(self):
        fd = gzip.open(self.fnm, 'rb')
        self.lines = fd.readlines()
        fd.close()

    def parseLines(self):
        for line in self.lines:
            line = str(line, 'UTF-8').strip()
            if line.startswith('@'):
                if ('read' in line) and ('barcode' in line):
                    self.hdr.append(line[1:])
                    self.stag = True
            else:
                if self.stag:
                    self.seq.append(line)
                    self.stag = False

    def Process(self):
        self.readLines()
        self.parseLines()

    def GetHdrs(self):
        return self.hdr

    def GetSeqs(self):
        return self.seq

    def GetNumRecords(self):
        return len(self.seq)
