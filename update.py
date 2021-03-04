import numpy as np
import sys
def obtainnatoms(filedftinput):
  f=open(filedftinput,'r');
  lines=f.readlines();
  for i in lines:
    if i.find("nat")!=-1:
      arra=i.split('=');
      arra=arra[1].split(',');
      natom=int(arra[0]);
  return natom;
def writenewscf(natoms,scfin,atomposition,filename):
    scffiles=open(scfin,'r');
    newfilename=open(filename,'w');
    lines=scffiles.readlines();
    tick1=0;
    tick2=0;
    skiplines=0;
    for i in range(len(lines)):
      if lines[i].find("ATOMIC_POSITIONS (angstrom)")!=-1:
        newfilename.write("ATOMIC_POSITIONS (angstrom)\n");
        tick2=i;
        for j in range(natoms):
          temp="";
          for t in range(3):
            temp=temp+" "+'{:12.8f}'.format(atomposition[j,t]);
          newfilename.write(lines[j+i+1].split()[0]+" "+temp+"\n");
        skiplines=natoms+skiplines+1;
      if skiplines > 0:
        skiplines=skiplines-1;
        continue;
      else:
        newfilename.write(lines[i]);
    newfilename.close();
    scffiles.close();
def updatedft(dftinput,dftoutput):
    dftainput='DFTAFTER.in';
    natoms=obtainnatoms(dftinput)
    filere=open(dftoutput,'r');
    lines=filere.readlines();
    Endsignal=0;
    atomsp=np.zeros((natoms,3))
    length=len(lines);
    for i in range(length):
      if lines[i].find("End final")!=-1:
        Endsignal=1;
      if lines[i].find("ATOMIC_POSITIONS")!=-1:
        for j in range(natoms):
          line=lines[i+j+1].split();
          for k in range(3):
            atomsp[j][k]=float(line[k+1]);
    writenewscf(natoms,dftinput,atomsp,dftainput)
updatedft(sys.argv[1],sys.argv[2])
