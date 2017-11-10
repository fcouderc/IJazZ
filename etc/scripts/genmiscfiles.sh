#!/bin/bash

# Quick and dirty generation of root_include_path.out and asap_config.h

LANG=C
CXXFLAGS=`cat`
echo "${CXXFLAGS}" |tr ' ' '\n' |grep '^-I' |sed -e 's/^-I//' \
  | while read i; do readlink -f $i; done |xargs |tr ' ' ':' \
  > tmp.out

if ! cmp --silent tmp.out root_include_path.out; then
  mv tmp.out root_include_path.out
fi

echo '#ifndef ASAP_CONFIG_H' > tmp.out
echo "${CXXFLAGS}" |tr ' ' '\n'|grep '^-D'|sort|uniq \
  |while read i; do
  def=`echo "$i" | sed -e 's/^-D//; s/=.*//'`
  if [[ "${i}" == *=* ]]; then
    target=`echo "${i}" | sed -e 's/.*=//'`
  else
    target=1
  fi
  echo "#ifndef ${def}"
  echo "#define ${def} ${target}"
  echo "#endif"
done >> tmp.out
echo "#endif" >> tmp.out

if ! cmp --silent tmp.out asap_config.h; then
  mv tmp.out asap_config.h
fi

rm -f tmp.out

