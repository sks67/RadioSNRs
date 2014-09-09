#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create the configuration file for input       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import ConfigParser

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ENTER SECTION AND NAMES HERE       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
section1 = 'ObsParameters'
section2 = 'InputFiles'
path_to_folder = '/Users/sumits2k/Desktop/Research/SNResearch2/RADSNRS/Inputs/'
hmap = 'LMC60.M0NHC3.FITS' 
snrlums = 'lmcradiolumscopy.txt'
snrsizes = 'lmc_angdiams.txt'
snrdens = 'lmc_snrs_opaccopy.txt'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Config = ConfigParser.ConfigParser()

#Create the config file
cfgfile = open('config.ini','w')

#Add Sections to the file

Config.add_section(section1)
Config.set(section1,'Distance',0.05)
Config.set(section1,'NumberofSNRs',54)
Config.add_section(section2)
Config.set(section2,'Path',path_to_folder)
Config.set(section2,'HydrogenMap',hmap)
Config.set(section2,'Luminosity',snrlums)
Config.set(section2,'Size',snrsizes)
Config.set(section2,'Density',snrdens)
Config.write(cfgfile)
cfgfile.close()
