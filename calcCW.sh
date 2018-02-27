export DISPLAY=:0
cd ~/Documents/MCCorradoExt

sleep 240

# FIND FILE NAME AND CREATE temp.R
(cat mcCW.R & head -n 1 logs/workCW.log) > temp.R

# Check if an R process is already running
isrunning=$((top -b -n 1 -o %CPU | grep ^\\s[0-9] | head -n 12 | grep tgauss$) && echo yes || echo no)


# If R is not running send to R and when work is done delete the first line from tempCW.log
[[ $isrunning == no ]] && R -f temp.R && ((tail -n +2 logs/workCW.log > logs/tempCW.log) | mv logs/tempCW.log logs/workCW.log)

if /usr/bin/zenity --question --timeout=10 --title="Calculation ended" --text="System is going to reboot. Do you want to abort it?"
then /usr/bin/zenity --info --text="Reboot is aborted"
else shutdown -r now 
fi
