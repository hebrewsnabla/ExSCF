hf=$1
for i in 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.7 1.8 2.0 2.2 2.5
do
  cp ${hf}_tpl.gjf ${hf}-${i}.gjf
  sed -i "s/dis/$i/g" ${hf}-${i}.gjf
  g09 ${hf}-${i}.gjf
done
