# export DISPLAY=:0
cd ~/Documents/MCCorradoExt

sleep 240

# FIND FILE NAME AND CREATE temp.R
(cat mcCWplus.R & head -n 1 logs/workCW.log) > temp.R

# Check if an R process is already running
isrunning=$((top -b -n 1 -o %CPU | grep ^\\s[0-9] | head -n 12 | grep tgauss$) && echo yes || echo no)


# If R is not running send to R and when work is done delete the first line from tempCW.log
[[ $isrunning == no ]] && R -f temp.R && ((tail -n +2 logs/workCW.log > logs/tempCW.log) | mv logs/tempCW.log logs/workCW.log)

sleep 120 && bash calcCW.sh
