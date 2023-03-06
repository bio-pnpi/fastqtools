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

class ExpandPrimer(object):
    def __init__(self, primer):
        self.primerSize = len(primer)
        self.forwardPrimers = [primer]
        self.reversePrimers = []
        self.degeneratedSubstitutions = {
            "M": ["A", "C"],
            "R": ["A", "G"],
            "W": ["A", "T"],
            "S": ["C", "G"],
            "Y": ["C", "T"],
            "K": ["G", "T"],
            "V": ["A", "C", "G"],
            "H": ["A", "C", "T"],
            "D": ["A", "G", "T"],
            "B": ["C", "G", "T"],
            "N": ["A", "C", "G", "T"]
        }
        self.complimentaryNucleotides = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }

    def ExpandPosition(self, position):
        primers = []
        for primer in self.forwardPrimers:
            if primer[position] in self.degeneratedSubstitutions.keys():
                for degPos in self.degeneratedSubstitutions[primer[position]]:
                    self.forwardPrimers.append(primer[:position] + degPos + primer[position + 1:])
            else:
                primers.append(primer)
        self.forwardPrimers = sorted(set(primers))

    def ExpandAllPositions(self):
        for i in range(self.primerSize):
            self.ExpandPosition(i)
        if self.forwardPrimers is not None:
            self.forwardPrimers = sorted(set(self.forwardPrimers))

    def ReversePrimers(self):
        for primer in self.forwardPrimers:
            for i in range(self.primerSize):
                primer = primer[:i] + self.complimentaryNucleotides[primer[i]] + primer[i + 1:]
            self.reversePrimers.append(primer[::-1])

    def Expand(self):
        self.ExpandAllPositions()
        self.ReversePrimers()

    def GetForward(self):
        return self.forwardPrimers

    def GetReverse(self):
        return self.reversePrimers
