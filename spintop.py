import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sb
from match import get_transform, quaternion_to_euler_angles,\
    find_orthonormals, quaternion_twist, quaternion_to_axis_angle


reportString = '''
% !TEX TS-program = pdflatex
% !TEX encoding = UTF-8 Unicode


\\documentclass[11pt]{{scrartcl}} 

\\usepackage[utf8]{{inputenc}} 

\\usepackage{{geometry}} 
\\geometry{{a4paper}} 

\\usepackage{{graphicx}} 
\\usepackage{{subfig}} 


\\usepackage{{sectsty}}
\\allsectionsfont{{\\sffamily\\mdseries\\upshape}} 

\\title{{{title}}}
\\subtitle{{Oscillation Range Report}}


\\begin{{document}}
\\setlength{{\\footskip}}{{3cm}}
\\maketitle

\\section{{Completeness and $I/\\sigma$}}

{completeness}
{iOverSigma}



\\newpage
\\section{{Oscillation Range Fits for Different Values}}
{fits}



\\end{{document}}

'''

class RangeReport(object):

    def __init__(self, dataName):
        self.comp = None
        self.iOverSigma = None
        self.fits = []
        self.dataName = dataName.replace('_', '\\_')
        self.startingValue = None

    def setCompletenessFigure(self, name):
        self.comp = name

    def setIOverSigmaFigure(self, name):
        self.iOverSigma = name

    def addFit(self, name):
        self.fits.append(name)

    def setBaseValue(self, value):
        self.startingValue = value

    def report(self):
        string = reportString.format(completeness='\\includegraphics[width=1\\textwidth]{{{}}}\n'.format(self.comp),
                                     iOverSigma='\\includegraphics[width=1\\textwidth]{{{}}}\n'.format(self.iOverSigma),
                                     fits='\n'.join(['\\includegraphics[width=.5\\textwidth]{{{{{}}}.eps}}'.format(fit[:-4]) for fit in self.fits]),
                                     title='{} at Oscillation Range {:6.4f}'.format(self.dataName, self.startingValue))
        with open('./plots/report.tex', 'w') as fp:
            fp.write(string)

class XDSIntegrateLog(object):
    def __init__(self, fileName, overrideTitle=None, savePlotsTo=None, fileNameCallback=None):
        self.fileName = fileName
        self.range = None
        self.imageCacheSize = None
        self.rotationAxis = None
        self.oscillationRange = None
        self.startingAngle = None
        self.numberOfBlocks = None
        self.blocks = []
        self.imagesInBlock = []
        self.overrideTitle = overrideTitle
        self.savePlotsTo = savePlotsTo
        self.filNameCallback = fileNameCallback

    def read(self):
        currentBlock = -1
        with open(self.fileName, 'r') as fp:
            for line in fp.readlines():
                line = line.strip()
                if line.startswith('!'):
                    continue
                if line.startswith('DATA_RANGE'):
                    self.range = tuple([int(word) for word in line.split()[-2:] if word])
                if line.startswith('NUMBER_OF_IMAGES_IN_CACHE'):
                    self.imageCacheSize = int(line.split()[-1])
                    self.numberOfBlocks = (self.range[1] - self.range[0]) // self.imageCacheSize + (
                        1 if (self.range[1] - self.range[0]) % self.imageCacheSize else 0)
                if line.startswith('ROTATION_AXIS'):
                    self.rotationAxis = np.array([float(word) for word in line.split()[-3:] if word])
                if line.startswith('OSCILLATION_RANGE'):
                    self.oscillationRange = float(line.split()[-2])
                if line.startswith('STARTING_ANGLE'):
                    self.startingAngle = float(line.split()[1])

                if 'PROCESSING OF IMAGES' in line:
                    line = [word for word in line.split() if word]
                    self.imagesInBlock.append(int(line[5]) - int(line[3]))

                    currentBlock += 1
                    if currentBlock > 0:
                        block = XDSIntegrateBlock(currentBlock-1, a, b, c, self.rotationAxis)
                        self.blocks.append(block)
                    a = None
                    b = None
                    c = None
                    continue

                if line.startswith('COORDINATES OF UNIT CELL A-AXIS'):
                    a = np.array([float(word) for word in line.split()[-3:] if word])
                if line.startswith('COORDINATES OF UNIT CELL B-AXIS'):
                    b = np.array([float(word) for word in line.split()[-3:] if word])
                if line.startswith('COORDINATES OF UNIT CELL C-AXIS'):
                    c = np.array([float(word) for word in line.split()[-3:] if word])
            block = XDSIntegrateBlock(currentBlock, a, b, c, self.rotationAxis)
            self.blocks.append(block)


        onorm, _ = find_orthonormals(self.rotationAxis)

        twists = []
        frameNumbers = []
        frameNumber = self.imagesInBlock[0]
        angles = []
        for i, b1 in enumerate(self.blocks[1:]):
            q = b1 - self.blocks[0]
            twist = quaternion_twist(q, self.rotationAxis, onorm)
            twists.append(twist)
            frameNumber += self.imagesInBlock[i+1]
            frameNumbers.append(frameNumber)
            angles.append(self.oscillationRange*frameNumber)

        A = np.vstack([angles, np.ones(len(angles))]).T
        m, b = np.linalg.lstsq(A, twists, rcond=None)[0]
        x = angles
        y = [i*m+b for i in x]
        sb.set_style("darkgrid")
        plt.plot(angles, twists, label='Offset About Rotation Axis')
        m = m * self.oscillationRange
        plt.plot(x, y, label='recommended correction = {:7.5f} deg/frame'.format(m))
        plt.plot(x, [t-i for t,i in zip(twists, y)], label='expected offset after correction')
        plt.legend(loc=2)
        plt.title('/'.join(self.fileName.split('/')[-5:-1]) if not self.overrideTitle else self.overrideTitle)
        plt.xlabel('Rotation Angle in degrees')
        plt.ylabel('Angular Offset in degrees')

        print('Old OSCILLATION_RANGE= {:9.7f}'.format(self.oscillationRange))
        print('Use OSCILLATION_RANGE= {:9.7f}'.format(self.oscillationRange+m))
        if self.savePlotsTo:
            plt.gcf()
            saveName = os.path.join(self.savePlotsTo ,'fit_{:6.4f}.eps'.format(self.oscillationRange))
            plt.savefig(saveName)
            self.filNameCallback('fit_{:6.4f}.eps'.format(self.oscillationRange))
        plt.show()





class XDSIntegrateBlock(object):
    def __init__(self, i, a, b, c, rotationAxis=None):
        self.i = i
        self.a = a
        self.b = b
        self.c = c
        self.points = (a, b, c)
        self.A = a/np.linalg.norm(a)
        self.B = b / np.linalg.norm(b)
        self.C = c / np.linalg.norm(c)
        self.Points = (self.A, self.B, self.C)


    def __sub__(self, other):
        t = get_transform(self.points, other.points, matrix=False)
        return t



def test(baseDir):
    import os

    for root, dir, files in os.walk(baseDir):
        # if 'bkp' in root or 'mods' in root:
        #     continue
        if 'bkp' in root:
            continue
        if 'brute' in root:
            continue
        # if not 'x05' in root:
        # if not 'sil1_boxE2_x15_42' in root:
        if not 'id13' in root:
            continue
        if 'INTEGRATE.LP' in files:
            print()
            print(root)
            try:
                fileName = os.path.join(root, 'INTEGRATE.LP')
                log = XDSIntegrateLog(fileName)
                log.read()
            except:
                print('Something went wrong.')
    # plt.show()


def bruteForce(baseDir, name):
    import os
    import numpy as np
    import subprocess

    saveDir = os.path.join(baseDir, 'plots')
    os.makedirs(saveDir)

    cwd = os.getcwd()
    os.chdir(baseDir)

    report = RangeReport(name)

    baseLine = None
    baseValue = None

    devnull = open(os.devnull, 'w')
    comps, iovers = [], []
    compsH, ioversH = [], []
    compsL, ioversL = [], []
    values = []
    for offset in np.arange(-0.002, 0.003, 0.0001, float):
        try:
            os.remove(os.path.join(baseDir, 'INTEGRATE.LP'))
        except FileNotFoundError:
            pass
        with open(os.path.join(baseDir, 'XDS.INP')) as fp:
            lines = fp.readlines()
            for i, line in enumerate(lines):
                if 'OSCILLATION_RANGE' in line:
                    if not baseLine:
                        baseLine = (i,line)
                    start, value = line.strip().split('=')
                    value = float(value.strip().split()[0])
                    if not baseValue:
                        baseValue = value
                        report.setBaseValue(baseValue)
                    value = baseValue + offset
                    print('\nTrying Oscillation Range of {:8.4f}'.format(value))
                    lines[i] = '{}={:8.6f} \n'.format(start, value)
                    # print(lines[i])
                    break
        with open(os.path.join(baseDir, 'XDS.INP'), 'w') as fp:
            fp.write(''.join(lines))

        subprocess.call(['xds_par'], stdout=devnull)
        fileName = os.path.join(baseDir, 'INTEGRATE.LP')
        log = XDSIntegrateLog(fileName, overrideTitle='OSCILLATION_RANGE = {:6.4f}'.format(baseValue + offset), savePlotsTo=saveDir, fileNameCallback=report.addFit)
        try:
            log.read()
        except FileNotFoundError:
            print('Not stable')
        except TypeError:
            print('Unknown Error')
        else:
            try:
                comp, Iover, allComps, allIovers = getCompletenessAndIoverSigma(baseDir)
            except:
                print('Cannot Parse CORRECT.LP')
            else:
                comps.append(comp)
                iovers.append(Iover)
                values.append(value)

                compsH.append(allComps.pop())
                compsL.append(allComps.pop(0))

                ioversH.append(allIovers.pop())
                ioversL.append(allIovers.pop(0))
    with open(os.path.join(baseDir, 'XDS.INP'), 'w') as fp:
        lines[baseLine[0]] = baseLine[1]
        fp.write(''.join(lines))

    plt.plot(values, comps, label='Overall Completeness')
    plt.plot(values, compsH, label='Completeness in highest resolution shell')
    plt.plot(values, compsL, label='Completeness in lowest resolution shell')
    plt.title('Completeness Versus Oscillation Range')
    plt.xlabel('Oscillation Range in degrees')
    plt.ylabel('Completeness')
    plt.legend(loc=3)
    plt.savefig(os.path.join(saveDir, 'completeness.eps'))
    plt.show()
    report.setCompletenessFigure('completeness.eps')

    plt.plot(values, iovers, label='Overall I/Sigma')
    plt.plot(values, ioversH, label='I/Sigma in highest resolution shell')
    plt.plot(values, ioversL, label='I/Sigma in lowest resolution shell')
    plt.title('I/Sigma Versus Oscillation Range')
    plt.xlabel('Oscillation Range in degrees')
    plt.ylabel('I/Sigma')
    plt.legend(loc=3)
    plt.savefig(os.path.join(saveDir, 'iOverSigma.eps'))
    plt.show()
    report.setIOverSigmaFigure('iOverSigma.eps')




    report.report()
    os.chdir(cwd)


def getCompletenessAndIoverSigma(baseDir):
    completeness = []
    IoverSigmas = []
    with open(os.path.join(baseDir, 'CORRECT.LP')) as fp:
        switch = False
        switch2 = False
        for i, line in enumerate(fp.readlines()):
            if switch and switch2 and i>switch2+1:# and 54 < i - switch < 65:
                line = [word for word in line.strip().split() if word]
                if not line:
                    break
                try:
                    complete, IoverSigma = float(line[4][:-1]), float(line[8])
                except ValueError:
                    break
                completeness.append(complete)
                IoverSigmas.append(IoverSigma)

            elif 'STATISTICS OF SAVED DATA SET' in line:
                switch = i
            elif switch and 'LIMIT     OBSERVED  UNIQUE  POSSIBLE' in line:
                switch2 = i
    comp = completeness.pop()
    IoverSigma = IoverSigmas.pop()
    return comp, IoverSigma, completeness, IoverSigmas


if __name__ == '__main__':
    test('/home/jens/tg/OscillationWidth/')
    # bruteForce('/home/jens/tg/OscillationWidth/brute/ir7b_pos14_x1/pos14_x1/id14', 'ir7b_id14')
    # print(getCompletenessAndIoverSigma('/home/jens/tg/OscillationWidth/brute/ir7b_pos14_x1/pos14_x1/id13'))
    # exit()
    # bruteForce('/home/jens/tg/OscillationWidth/brute/ir7b_pos14_x1/pos14_x1/id13', 'ir7b_id13')
    # bruteForce('/home/jens/tg/OscillationWidth/brute/sil1_x05_17/x05_17', 'sil1_x05_17')
    # test('/home/jens/tg/data/bkp/data/')
