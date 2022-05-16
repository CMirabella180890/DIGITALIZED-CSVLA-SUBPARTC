# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 16:04:18 2022

@author: claum
"""
import os

# GETTING MAIN PATH
main_path = os.getcwd()
print(" Working directory:\n ", main_path, "\n")
####################################
##### CHANGE WORKING DIRECTORY #####
####################################
init_dir = main_path + '\initialization'
os.chdir(init_dir)
print(" Working directory:\n ", init_dir, "\n")
################################
##### INITIALIZATION PHASE #####
################################
print(" +++ INITIALIZATION PHASE +++ \n ")
exec(open("initialization.py").read())
####################################
##### CHANGE WORKING DIRECTORY #####
####################################
print(" Returning to main directory.\n")
os.chdir(main_path)
print(" Working directory:\n ", main_path)
####################################
##### CHANGE WORKING DIRECTORY #####
####################################
csvla_dir = main_path + '\csvla'
os.chdir(csvla_dir)
print(" Working directory:\n ", csvla_dir, "\n")
################################
##### MANOEUVRING ENVELOPE #####
################################
print(" +++ FLIGHT ENVELOPE DIAGRAM +++ \n ")
exec(open("CalcFlightEnvelope.py").read())
exec(open("CalcGustEnvelope.py").read())
exec(open("CalcFinalEnvelope.py").read())
exec(open("CalcBalancingLoads.py").read())
exec(open("CalcInterpDistrLiftDragPitchingMom.py").read())
exec(open("CalcShearBendingTorsion.py").read())