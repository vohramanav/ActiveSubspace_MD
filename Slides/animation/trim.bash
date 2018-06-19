#!/bin/bash

p = 1
for i in {0..1280..20}
  do
    if [ $i -lt 10 ]
    then
      convert snap000$i.png -trim snap000$i.png
      mv snap000$i.png snap$p.png
      let p++
    fi
    
    if [ $i -ge 10 ] && [ $i -lt 100 ]
    then
      convert snap00$i.png -trim snap00$i.png
      mv snap00$i.png snap$p.png
      let p++
    fi

    if [ $i -ge 100 ] && [ $i -lt 1000 ]
    then
      convert snap0$i.png -trim snap0$i.png
      mv snap0$i.png snap$p.png
      let p++
    fi

     if [ $i -ge 1000 ] && [ $i -lt 10000 ]
    then
      convert snap$i.png -trim snap$i.png
      mv snap$i.png snap$p.png
      let p++
    fi

  done

  convert -loop 10 *.png md_si.gif
  



