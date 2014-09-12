#
# This is a common set of routines that IPython scripts accessing the database use.
#


# Connect to the database - currently on my office computer. If not running this directly on that computer,
#  an SSH tunnel can be set up via PuTTy or BitviseSSH (along with other tools) using the standard port 3306.
#import MySQLdb as mdb
import pymysql as mdb
import numpy as np
import sys

try:
	db = mdb.connect(host="localhost", user="root", passwd="koW4Gjgot2", db="psh_data")
except mdb.Error, e:
	sys.stderr.write("[ERROR] %d: %s\n" % (e.args[0], e.args[1]))
	exit()
#cursor = db.cursor(mdb.cursors.DictCursor)
cursor = db.cursor()


def GetKappaTablenames(cursor, mysqlstring, points, minK=0.0, maxK=1000000.0):
    """ The psh_data database contains entries detailing attributes of runs with different kappa values.
        This retrieves the names of all of the data tables. minK and maxK allow for selection of
        a smaller range than the full set. """
    cursor.execute("SELECT kappa,datatable FROM phases WHERE %s AND quadpoints='%s'" % (mysqlstring, points))
    #numrows = int(cursor.rowcount)
    entries = []
    rows = cursor.fetchall()
    for row in rows:
        if row[0] < minK or row[0] > maxK:
            continue
        entries.append(row)
    entries.sort()  # Sorts in increasing order of kappa
    return entries


def GetSLS(cursor, mysqlstring, lvalue, points, minK=0.0, maxK=1000000.0):
    """ Retrieves the SLS entries for all data tables. minK and maxK allow for selection of
        a smaller range than the full set. """
    #@TODO: Combine this with GetKappaTablenames?
    slstable = []
    kappatable = []
    
    cursor.execute("SELECT kappa,sls FROM phases WHERE %s AND quadpoints='%s'" % (mysqlstring, points))
    entries = []
    rows = cursor.fetchall()
    for row in rows:
        if row[0] < minK or row[0] > maxK:
            continue
        entries.append(row)
    entries.sort()  # Sorts in increasing order of kappa
    
    for entry in entries:
        kappatable.append(entry[0])
        slstable.append(entry[1])
    kappa = np.array(kappatable)
    phases = -np.array(slstable)
    cross = 4.0 * (2*lvalue + 1) * np.sin(phases)**2 / kappa**2
    
    return [phases, cross]


def GetPhaseShifts(cursor, entries, numterm, lvalue, kohnmethod):
    """ This retrieves the phase shift for each data table selected at N (numterm), for the Kohn method given
        (e.g. S-matrix). This constructs a kappa list as well, showing what kappa values are present. These
        are also converted to energies and cross sections. """
    phase = []
    kappa = []
    for i in range(len(entries)):
        tablename = entries[i][1]  # Select the phase shift at N terms
        cursor.execute("SELECT %s FROM %s WHERE n = %i" % (kohnmethod, tablename, numterm))
        s = cursor.fetchone()  # Should *only* ever have one entry
        if s == None:  # Some may be incomplete - do not attempt to use them
            print tablename + " does not have an n = %i entry." % numterm
            continue
        kappa.append(entries[i][0])
        phase.append(s[0])
    phase = np.array(phase)
    kappa = np.array(kappa)
    energy = 27.21138505/4.0 * kappa**2
    cross = 4.0 * (2*lvalue + 1) * np.sin(phase)**2 / kappa**2

    return [kappa, energy, phase, cross]


def GetPhaseArray(cursor, entries, numterms, lvalue, kohnmethod):
    """ Similar to GetPhaseShifts but does this for all N values """
    phase = []
    kappa = []
    cross = []

    for i in range(len(entries)):
        tablename = entries[i][1]
        cursor.execute("SELECT %s FROM %s" % (kohnmethod, tablename))
        p = []
        rows = cursor.fetchall()
        if len(rows) < numterms+1:
            print tablename + " does not have the correct number of entries."
            continue
        #for row in rows:
        for j in range(numterms+1):  # May want less than all possible terms
            row = rows[j]
            p.append(row[0])
        k = entries[i][0]
        kappa.append(entries[i][0])
        phase.append(p)
        cross.append(4.0 * (2*lvalue + 1) * np.sin(np.array(p))**2 / k**2)

    kappa = np.array(kappa)
    phase = np.array(phase).T  # Rearrange to sets of N
    cross = np.array(cross).T
    energy = 27.21138505/4.0 * kappa**2
       
    return [kappa, energy, phase, cross]


def GetSubset(kappafull, kappapoints, phasecross):
    """ Creates a subset of the phase shift or cross section lists, using only the kappa values in kappapoints """
    k = kappafull.tolist()
    kindex = []
    for i,j in enumerate(kappapoints):
        kindex.append(k.index(j))
        
    return phasecross[kindex]
