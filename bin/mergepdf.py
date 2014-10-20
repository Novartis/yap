#!/usr/bin/env python
# Copyright (C) 2010 Milinda Pathirage

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
from pyPdf import PdfFileWriter, PdfFileReader
import glob


def mergePDFFiles(outputFile, filesToBeMerged):
    output = PdfFileWriter()

    if(len(filesToBeMerged) == 0):
        print 'Empty Input File List'
        return

    for inFile in filesToBeMerged:
        # print 'Adding file' + inFile + ' to the out put'
        # Read the input PDF file
        input = PdfFileReader(file(inFile, "rb"))
        # Add every page in input PDF file to output
        for page in input.pages:
            output.addPage(page)
    # print 'Writing the final out put to file system'
    # Out put stream for output file
    outputStream = file(outputFile, "wb")
    output.write(outputStream)
    outputStream.close()


if __name__ == '__main__':
    i = 0
    outputFile = 'merged.pdf'

    # CR559066_1
    #inputFiles = glob.glob('CR559066_1/no*/post*/*.txt')
    #inputFiles = ["junk.f90.pdf","junk.f90.pdf"]
    for arg in sys.argv:
        i = i + 1

        # Getting out file
        if arg == '-o':
            outputFile = sys.argv[i]
            # print 'Output File: ' + outputFile

        # Extracting Input files
        if arg == '-i':
            outfileOptionPos = sys.argv.index('-o')
            if i < outfileOptionPos:
                inputFiles = sys.argv[i: outfileOptionPos]
                filesStr = ",".join(inputFiles).replace(",", " ")
                # print 'Input Files: ' + filesStr
            else:
                inputFiles = sys.argv[i:]
                filesStr = ",".join(inputFiles).replace(",", " ")
                # print 'Input Files: ' + filesStr

    # Merging PDF files
    # print 'Merging PDF Files......'
    mergePDFFiles(outputFile, inputFiles)
