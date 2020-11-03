hf=$1
g09m="g09 -exedir=/share/home/srwang/g09m/exe-dir:$GAUSS_EXEDIR"
for i in 0.75 0.80  0.90 0.95 1.00 1.05 1.10 1.15 1.20 # 1.30 1.40 1.50 1.70 2.00 2.20 2.50
do
  cp ${hf}_tpl.gjf ${hf}-${i}.gjf
  sed -i "s/dis/$i/g" ${hf}-${i}.gjf
  $g09m ${hf}-${i}.gjf
done
