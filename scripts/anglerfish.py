import subprocess
import sys

print("!!!\nanglerfish.py is deprecated. Use the command 'anglerfish' from version 0.4.2 on \n!!!")
args = ' '.join(sys.argv[1:])
p = subprocess.Popen(f"anglerfish {args}", stdout=subprocess.PIPE, shell=True)
out, err = p.communicate() 
result = out.split(b'\n')
for lin in result:
        print(str(lin, 'UTF-8'))