#!/usr/bin/python
#-----------------------------------------------------------------------------------------------#
# Name - metanalysis.py										#
# Desc - Untargeted metabolomics annotation.							#
# Author - Difei Wang (dw670@georgetown.edu); Lin An (la512@georgetown.edu)			#
#-----------------------------------------------------------------------------------------------#

# Standard dependancies
from optparse   import OptionParser
from sys        import stderr
from subprocess import Popen
from subprocess import PIPE
import re
import csv
import MySQLdb
import operator

db  = MySQLdb.connect(host  = "localhost", # your host, usually localhost
                     user   = "root", # your username
                     passwd = "", # your password
                     db     = "chempro") # name of the database
# Catch user input options
def parseArguments ():
    parser = OptionParser ( "Usage: %prog -f <m/zlist> -p <positive> -n <negative>" )
    parser.add_option ( "-f", "--file", dest = "Monoisotopic_weight", help = "The report file of differential expressed chemicals." )
    parser.add_option ( "-m", "--ppm", type="int", dest = "PPM", help = "The ppm value" )
    parser.add_option ( "-p", "--positive", action="store_true", dest = "Positive", help = "The positive polarity" )
    parser.add_option ( "-n", "--negative", action="store_true", dest = "Negative", help = "The negative polarity" )
    parser.add_option ( "-i", "--isoform", action="store_true", dest = "Isoform", help = "Find Isoforms" )
    parser.add_option ( "-a", "--adduct", action="store_true", dest = "Adduct", help = "Find Adduct" )
    (options, arguments) = parser.parse_args ()
    if ( not options.Monoisotopic_weight ):
        parser.error ( "Requires a m/z value list of chemical. Run with --help for help." )
        quit ( 1 )
    if ( not options.PPM ):
        parser.error ( "Requires a PPM value for database search. Run with --help for help." )
        quit ( 1 )
    if ( not options.Positive ) and ( not options.Negative ):
        parser.error ( "Requires a positive definition or a negative definition for analysis. Run with --help for help." )
    if ( options.Positive ) and ( options.Negative ):
        parser.error ( "Requires just one definition for analysis. Run with --help for help." )
    return options

# Parse user input files
def parsefiles (inputfile, isot):
    cfile = open(inputfile)
    data  = csv.DictReader(cfile)
    querylist = []
    for r in data:
        iso = r["isotopes"]
        if isot:
            ma = re.match(r'\[\d+\]\[M([\+\d]*)\]([\d]*)\+',iso)
            if ma:
                if ma.group(1) == "":
                    m1 = 0
                else:
                    m1 = float(ma.group(1)[1])
                if ma.group(2) == "":
                    m2 = 1
                else:
                    m2 = float(ma.group(2))
                mass = str(float(r["mzmed"])*m2 - m1)
                querylist.append((mass,iso, r["fold"], r["tstat"], r["pvalue"]))
            mb = re.match(r'\[\d+\]\[M\]\+',iso)
            if mb:
                continue
            else:
                querylist.append((r["mzmed"], "[M]+", r["fold"], r["tstat"], r["pvalue"]))
        else:
            querylist.append((r["mzmed"], "[M]+", r["fold"], r["tstat"], r["pvalue"]))
    cfile.close()
    return querylist

# Adducts list
def adductlist (init=False):
    if init:
        data = [("[M+H]+", 1.0073),
                ("[M+Na]+", 22.98977),
                ("[M+K]+", 38.961708)]
    else:
        data = [("[M-H]-", -1.0073)]
    return data

# Print Result
def printdb (cur,querylist,ppm,adduct):
    mz     = float(querylist[0])
    interval = (mz)/(1.00/ppm*1000000)
    mass   = mz - adduct[1]
    top    = mass + interval
    bottom = mass - interval
    cur.execute("SELECT * FROM hmdbchem1 WHERE monisotopweight <= %s and monisotopweight >= %s",(top,bottom))
    for line in cur.fetchall():
        resultlist = [str(n).strip() for n in list(line)]
        realmass = float(resultlist[2])
        deltappm = abs(mz-realmass-adduct[1])/mz*1000000
        deltappm = round(deltappm,1)
        singleline = str(mz)+'\t'+adduct[0]+'\t'+querylist[1]+'\t'+'\t'.join(resultlist)+'\t'+str(deltappm)+'\t'+querylist[2]+'\t'+querylist[3]+'\t'+querylist[4]
        print singleline

# Run Query
if __name__ == "__main__":
    arguments = parseArguments ()
    cur = db.cursor()
    querylist = parsefiles( arguments.Monoisotopic_weight, arguments.Isoform )
    ppm = arguments.PPM
    if arguments.Positive:
        add = adductlist(init=True)
    if arguments.Negative:
        add = adductlist(init=False)
    for mz in querylist:
        for adduct in add:
            printdb(cur, mz, ppm, adduct)
