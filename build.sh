#!/bin/bash

javac -J-Xmx256m *.java
jar cf telomere.jar *.class
rm *.class

g++ find_telomere.c -o find_telomere -lz
