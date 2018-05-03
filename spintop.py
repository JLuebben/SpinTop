import os
import argparse
try:
    from matplotlib import pyplot as plt
except ImportError:
    PLOTAVAILABLE = False
else:
    PLOTAVAILABLE = True
# try:
#     import seaborn as sb
# except ImportError:
#     PLOTAVAILABLE = False


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

\\section{{$I/\\sigma$}}

{iover}



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

    def report(self, path):
        string = reportString.format(iover='\\includegraphics[width=1\\textwidth]{{{}}}\n'.format(self.iOverSigma),
                                     fits='\n'.join(['\\includegraphics[width=.5\\textwidth]{{{{{}}}.eps}}'.format(fit[:-4]) for fit in self.fits]),
                                     title='{} at Oscillation Range {:6.4f}'.format(self.dataName, self.startingValue))
        with open(os.path.join(path, 'report.tex'), 'w') as fp:
            fp.write(string)

class XDSIntegrateLog(object):
    def __init__(self, fileName, overrideTitle=None, savePlotsTo=None, fileNameCallback=None, savePlotsAs=None):
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
        self.savePlotsAs = savePlotsAs
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

    def fitOscillationCorrection(self):
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
        m = m * self.oscillationRange
        correctedTwists = [t-i for t,i in zip(twists, y)]

        self.plotData = {'angles': angles,
                         'twists': twists,
                         'correctionFunction': y,
                         'correctedTwists': correctedTwists,
                         'correction': m}
        return m

    def plot(self, show=False):
        if not self.plotData:
            raise ValueError('Plot data is None. Call "fitOscillationCorrection()" before calling "plot".')
        angles, twists, correctionFunction, correctedTwists = self.plotData['angles'], self.plotData['twists'], self.plotData['correctionFunction'], self.plotData['correctedTwists']
        m = self.plotData['correction']
        # sb.set_style("darkgrid")
        plt.plot(angles, twists, label='Offset About Rotation Axis')
        plt.plot(angles, correctionFunction, label='recommended correction = {:+6.4f} deg/frame --> {:6.4f}'.format(m, self.oscillationRange+m))
        plt.plot(angles, correctedTwists, label='expected offset after correction')
        plt.legend(loc=2)
        plt.title('/'.join(self.fileName.split('/')[-5:-1]) if not self.overrideTitle else self.overrideTitle)
        plt.xlabel('Rotation Angle in degrees')
        plt.ylabel('Angular Offset in degrees')

        # print('Old OSCILLATION_RANGE= {:6.4f}'.format(self.oscillationRange))
        # print('Use OSCILLATION_RANGE= {:6.4f}'.format(self.oscillationRange+m))
        if self.savePlotsTo:
            plt.gcf()
            saveName = os.path.join(self.savePlotsTo ,'fit_{:6.4f}.eps'.format(self.oscillationRange))
            plt.savefig(saveName)
            self.filNameCallback('fit_{:6.4f}.eps'.format(self.oscillationRange))

        if self.savePlotsAs:
            plt.savefig(self.savePlotsAs)
        if show:
            plt.show()
        plt.clf()

    def writeRawData(self, fileName):
        if not self.plotData:
            raise ValueError('Plot data is None. Call "fitOscillationCorrection()" before calling "writeRawData".')
        with open(fileName, 'w') as fp:
            fp.write('Angle, Twist, CorrectedTwist')
            for angle, twist, corrected in zip(self.plotData['angles'], self.plotData['twists'], self.plotData['correctedTwists']):
                fp.write('{},{},{}'.format(angle, twist, corrected))





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



def getCorrection(baseDir):
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


def fullScan(baseDir, name, plotDirName):
    import os
    import numpy as np
    import subprocess

    saveDir = os.path.join(baseDir, plotDirName)
    try:
        os.makedirs(saveDir)
    except FileExistsError:
        print('Directory <{}> already exists. Please use a different name.'.format(plotDirName))
        exit(2)

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
            log.fitOscillationCorrection()
            log.plot(show=False)
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

    # plt.plot(values, comps, label='Overall Completeness')
    # plt.plot(values, compsH, label='Completeness in highest resolution shell')
    # plt.plot(values, compsL, label='Completeness in lowest resolution shell')
    # plt.title('Completeness Versus Oscillation Range')
    # plt.xlabel('Oscillation Range in degrees')
    # plt.ylabel('Completeness')
    # plt.legend(loc=3)
    # plt.savefig(os.path.join(saveDir, 'completeness.eps'))
    # plt.show()
    # report.setCompletenessFigure('completeness.eps')

    plt.plot(values, iovers, label='Overall I/Sigma')
    plt.plot(values, ioversH, label='I/Sigma in highest resolution shell')
    plt.plot(values, ioversL, label='I/Sigma in lowest resolution shell')
    plt.title('I/Sigma Versus Oscillation Range')
    plt.xlabel('Oscillation Range in degrees')
    plt.ylabel('I/Sigma')
    plt.legend(loc=3)
    plt.savefig(os.path.join(saveDir, 'iOverSigma.eps'))
    # plt.show()
    report.setIOverSigmaFigure('iOverSigma.eps')




    report.report(plotDirName)
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


import numpy as np

def get_quaternion(lst1, lst2, matchlist):
    """
    Returns the quaternion representing the best possible
    transformation to minimize the distance between all
    points in 'cloud1' and 'cloud2'.
    The second return value is the fitting criteria
    representing how well both clouds match.
    (the larger the value the better.)
    """
    M = np.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    for i, coord1 in enumerate(lst1):
        x = np.matrix(np.outer(coord1, lst2[matchlist[i]]))
        M = M + x

    N11 = float(M[0][:, 0] + M[1][:, 1] + M[2][:, 2])
    N22 = float(M[0][:, 0] - M[1][:, 1] - M[2][:, 2])
    N33 = float(-M[0][:, 0] + M[1][:, 1] - M[2][:, 2])
    N44 = float(-M[0][:, 0] - M[1][:, 1] + M[2][:, 2])
    N12 = float(M[1][:, 2] - M[2][:, 1])
    N13 = float(M[2][:, 0] - M[0][:, 2])
    N14 = float(M[0][:, 1] - M[1][:, 0])
    N21 = float(N12)
    N23 = float(M[0][:, 1] + M[1][:, 0])
    N24 = float(M[2][:, 0] + M[0][:, 2])
    N31 = float(N13)
    N32 = float(N23)
    N34 = float(M[1][:, 2] + M[2][:, 1])
    N41 = float(N14)
    N42 = float(N24)
    N43 = float(N34)

    N = np.matrix([[N11, N12, N13, N14],
                   [N21, N22, N23, N24],
                   [N31, N32, N33, N34],
                   [N41, N42, N43, N44]])

    values, vectors = np.linalg.eig(N)
    w = list(values)
    mw = max(w)
    quat = vectors[:, w.index(mw)]
    quat = np.array(quat).reshape(-1, ).tolist()
    return quat, mw

def rotate_by_quaternion(cloud, quat):
    """
    Rotations the points in 'cloud' with the
    transformation represented by the quaternion
    'quat'.
    """
    rotmat = get_rotation_matrix_from_quaternion(quat)
    return rotate_by_matrix(list(cloud), rotmat)

def rotate_by_matrix(cloud, rotmat):
    """
    Rotates the points in 'cloud' by the transformation represented
    by the rotation matrix 'rotmat'.
    """
    return [np.array(np.dot(coord, rotmat).tolist())[0] for coord in cloud]

def get_rotation_matrix_from_quaternion(q):
    """
    Returns the rotation matrix equivalent of the given quaternion.

    This function is used by the get_refined_rotation() function.
    """
    R = np.matrix([[q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3],
                    2 * (q[1] * q[2] - q[0] * q[3]),
                    2 * (q[1] * q[3] + q[0] * q[2])],
                   [2 * (q[2] * q[1] + q[0] * q[3]),
                    q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3],
                    2 * (q[2] * q[3] - q[0] * q[1])],
                   [2 * (q[3] * q[1] - q[0] * q[2]),
                    2 * (q[3] * q[2] + q[0] * q[1]),
                    q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]]])
    return R

def get_transform(points1, points2, matchlist=None, use=3, matrix=False):
    """
    Returns the Quaternion/Matrix representation of
    of the transformation that transforms the points in
    'points1' to the coordinate system of 'points2'.

    The first three points in each list are used by default.
    It is assumed that the first point in list 1 matches
    the first point in list 2 and so forth...
    This behavior can be modified by using the 'matchlist'
    option that must be a list of integers specifying
    which points to map. 'matchlist=[2,1,0]' implies that
    the first element of 'points1' is mapped to the third
    element of 'points2' and so forth.

    The optional argument 'use' determines how many elements
    of the point lists should be used to determine the
    transformation. The default is three. A value of '0'
    implies the use of all elements.

    The boolean 'matrix' determines whether the
    transformation is returned in quaternion
    representation or in matrix representation.
    """
    if use == 0:
        use = len(points1)
    lst1 = points1[:use]
    lst2 = points2[:use]
    if not matchlist:
        matchlist = range(use)
    quat, _ = get_quaternion(lst1, lst2, matchlist)
    if not matrix:
        return quat
    else:
        return get_rotation_matrix_from_quaternion(quat)



import math
def quaternion_to_euler_angles(q):
    q0, q1, q2, q3 = q
    phi = np.arctan((2*(q0*q1 + q2*q3))/(1-2*(q1*q1 + q2*q2)))
    theta = np.arcsin(2*(q0*q2 - q3*q1))
    psi = np.arctan((2*(q0*q3 + q1*q2)) / (1-2*(q2*q2 + q3*q3)))

    phi = math.degrees(phi)
    theta = math.degrees(theta)
    psi = math.degrees(psi)
    return (phi, theta, psi)






rx90 = np.matrix(((1, 0, 0),
                  (0, 0, -1),
                  (0, 1, 0)))

ry90 = np.matrix(((0, 0, 1),
                  (0, 1, 0),
                  (-1, 0, 0)))

def find_orthonormals(vector):
    normal = vector/np.linalg.norm(vector)
    w = np.dot(normal, rx90).flatten().tolist()[0]
    d = np.dot(normal, w)
    if abs(d) > .6:
        w = np.dot(normal, ry90).flatten().tolist()[0]
    w = w/np.linalg.norm(w)
    onorm1 = np.cross(normal, w)
    onorm1 = onorm1 / np.linalg.norm(onorm1)

    onorm2 = np.cross(normal, onorm1)
    onorm2 = onorm2 / np.linalg.norm(onorm2)
    return onorm1, onorm2

def quaternion_twist(quat, axis, onorm):
    rotmat = get_rotation_matrix_from_quaternion(quat)
    transformed = np.array(np.dot(onorm, rotmat).tolist())[0]

    flattened = transformed - (np.dot(transformed, axis) * axis)
    flattened = flattened / np.linalg.norm(flattened)

    handedness = np.cross(onorm, flattened)
    # handedness = np.cross( flattened,onorm)
    handedness = handedness/np.linalg.norm(handedness)
    handedness = np.dot(handedness, axis)
    return math.degrees(np.arccos(np.dot(onorm, flattened)))*handedness


def quaternion_to_axis_angle(quat):
    v = quat[1:]
    return v/np.linalg.norm(v), math.degrees(2 *np.arccos(quat[0]))


def main():
    parser = argparse.ArgumentParser(description='Determine empirical correction factor for OSCILLATION_WITH parameter '
                                                 'in XDS integrations.')
    parser.add_argument('--path', '-p', type=str, default='./',
                        help="path to the directory containing the XDS.INP/INTEGRATE.LP/CORRECT.LP files.")
    parser.add_argument('--plot', '-P', action='store_true',
                        help='show a plot of the rotational offset.')
    parser.add_argument('--save', '-s', type=str, default='',
                        help="write plot to disk using the file name <SAVE>. Image file format is inferred from file"
                             " name.",)
    parser.add_argument('--scan', type=str, default='',
                        help='scan OSCILLATION_RANGE values around the the reported value by continuously calling '
                             'xds_par. This option might be useful if the data is particularly noisy or very finely '
                             'sliced. A directory named <SCAN> is created containing all generated plots as well as '
                             'a report.tex file that can be used to generate a diagnostic report via Latex.')
    parser.add_argument('--title', '-t', type=str, default='XDS_DataSet',
                        help='title used in the diagnostic report generated if the scan option is selected.')
    parser.add_argument('--raw', '-r', type=str, default='plotData.txt',
                        help='write the CSV file <RAW> containing the data to reproduce the plots.')
    args = parser.parse_args()
    if args.scan:
        fullScan(args.path, args.title, args.scan)
        return
    log = XDSIntegrateLog(os.path.join(args.path, 'INTEGRATE.LP'), savePlotsAs=args.save)
    log.read()
    m = log.fitOscillationCorrection()
    print('Old OSCILLATION_RANGE= {:6.4f}'.format(log.oscillationRange))
    print('Use OSCILLATION_RANGE= {:6.4f}'.format(log.oscillationRange + m))
    if args.plot or args.save:
        log.plot(show=args.plot)



if __name__ == '__main__':
    main()
    exit()
    # test('/home/jens/tg/OscillationWidth/')
    # bruteForce('/home/jens/tg/OscillationWidth/ir7b_pos14_x1/pos14_x1/id13', 'ir7b_id13')
    # bruteForce('/home/jens/tg/OscillationWidth/ir7b_pos14_x1/pos14_x1/id14', 'ir7b_id14')
    # bruteForce('/home/jens/tg/OscillationWidth/sil1_boxE2_x15_42/boxE2/x15_42', 'sil1_x15_42')
    # bruteForce('/home/jens/tg/OscillationWidth/brute/ir7b_pos14_x1/pos14_x1/id14', 'ir7b_id14')
    # print(getCompletenessAndIoverSigma('/home/jens/tg/OscillationWidth/brute/ir7b
    # _pos14_x1/pos14_x1/id13'))
    # exit()
    # bruteForce('/home/jens/tg/OscillationWidth/brute/ir7b_pos14_x1/pos14_x1/id13', 'ir7b_id13')
    fullScan('/home/jens/tg/OscillationWidth/brute/sil1_x05_17/x05_17', 'sil1_x05_17')
    # test('/home/jens/tg/data/bkp/data/')
