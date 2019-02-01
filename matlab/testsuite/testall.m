%% MAGCODE TEST SUITE


%% Start

clear all
clc

oldpath = path;
path(oldpath,'..')

%% Cuboids and other basics

test001a
test001b
test001c
test001d
test001e

testgrades01

test002a
test002b
test002c
test002d

test003a

%% Cylinders

testcyl01
testcyl02

%% Cuboid torques

testcuboidtorque01
testcuboidtorque02
testcuboidtorque03
testcuboidtorque04

%% End

path(oldpath)
