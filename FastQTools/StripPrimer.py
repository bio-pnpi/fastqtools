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

class StripPrimer(object):
    def __init__(self, seq, quality, fprimers, rprimers):
        self.forwardPrimers = fprimers
        self.reversePrimers = rprimers
        self.sequence = seq
        self.quality = quality
        self.minPrimerLength = 10
        self.stripped = []

    def StripForwardPrimerFromSeq(self, seq, quality):
        for fp in self.forwardPrimers:
            for fi in range(len(fp) - self.minPrimerLength):
                if seq.startswith(fp[fi:]):
                    if quality is not None:
                        return seq[len(fp[fi:]):], quality[len(fp[fi:]):]
                    else:
                        return seq[len(fp[fi:]):]
        if quality is not None:
            return seq, quality
        else:
            return seq

    def StripReversePrimerFromSeq(self, seq, quality):
        for rp in self.reversePrimers:
            for ri in range(self.minPrimerLength, len(rp)):
                if seq.endswith(rp[:ri]):
                    if quality is not None:
                        return seq[:len(seq) - len(rp[:ri])], quality[:len(seq) - len(rp[:ri])]
                    else:
                        return seq[:len(seq) - len(rp[:ri])]
        if quality is not None:
            return seq, quality
        else:
            return seq

    def Strip(self):
        if type(self.sequence) is list:
            for si in range(len(self.sequence)):
                self.sequence[si] = self.StripForwardPrimerFromSeq(self.sequence[si], None)
                self.sequence[si] = self.StripReversePrimerFromSeq(self.sequence[si], None)
        else:
             self.sequence, self.quality = self.StripForwardPrimerFromSeq(self.sequence, self.quality)
             self.sequence, self.quality = self.StripReversePrimerFromSeq(self.sequence, self.quality)

    def GetStripped(self):
        if self.quality is not None:
            return self.sequence, self.quality
        else:
            return self.sequence

