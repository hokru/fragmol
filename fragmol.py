import sys
import numpy as np
import argparse

"""

auto-fragmentation of non-cov bonded molecules
and by removing a bond

H.Kruse mail2holger@gmail.com 

"""
maxfrag=100
nat=0
elem = []
xyz = []
SYM=[]  # matrix to do symmetry operations

parser = argparse.ArgumentParser(description="make fragments of molecular clusters",epilog="42",usage='%(prog)s [options] <coordinate file>')

parser.add_argument("-cut", help="specify two atom numbers to cut a bond",type=int,nargs=2,metavar=("atom1","atom2"),default=[-1,-1])
parser.add_argument("molecule", help="molecular coordinate file (xyz format)",type=str,metavar="<coordinate file>")
args = parser.parse_args()

print('file              : ',args.molecule)
if args.debug:
   print("debugging mode turned on")




#! atomic radii from Mantina, Valero, Cramer, Truhlar "Atomic radii of elements"
# not used
rcov={'h':0.320,'he':0.370, 
    'li':1.300, 'be':0.990, 'b': 0.840, 'c': 0.750, 'n': 0.710,'o': 0.640,'f': 0.600,'ne':0.620,\
#    ! Na    Mg     Al     Si     P      S       Cl     Ar
    'na':1.60,'mg':1.400,'al':1.240,'si':1.140,'p':1.090,'s':1.040,'cl':1.000,'ar':1.010, \
#    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu  
    'k':2.000,'ca':1.740,'sc':1.590,'ti':1.480,'v':1.440,'cr':1.300,'mn':1.290,'fe':1.240,'co':1.180,'ni':1.170,'cu':1.220,
#    !  Zn     Ga     Ge     As     Se     Br    Kr
    'zn':1.200,'ga':1.230,'ge':1.200,'as':1.200,'se':1.180,'br':1.170,'kr':1.240, \
#    !  Rb    Sr
    'rb':2.150,'sr':1.900,
#    ! Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd
    'y':1.780,'zr':1.640,'nb':1.560,'mo':1.460,'tc':1.380,'ru':1.360,'rh':1.340,'pd':1.300,'ag':1.360,'cd':1.400,\
#     In    Sn      Sb      Te     I     Xe
    'in':1.420,'sn':1.400,'sb':1.400,'te':1.370,'i':1.320,'xe':1.360, \
#    ! Cs Ba
    'cs':2.380,'ba':2.060,\
#    ! La-Lu
#     1.94d0,1.84d0,1.90d0,1.73d0,1.86d0,1.85d0,1.83d0,1.82d0,1.81d0,1.80d0,1.79d0,1.77d0,1.77d0,1.78d0,1.74d0,  &
#    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    'hf':1.640,'ta':1.580,'w':1.500,'re':1.410,'os':1.360,'ir':1.320,'pt':1.300,'au':1.640,'hg':1.880,'ti':1.480,'pb':1.450,'bi':1.500,'po':1.420,'at':1.470,'rn':1.460}
#    ! Fr-Pu
#    2.42d0,2.11d0,2.01d0,1.90d0,1.84d0,1.83d0,1.80d0,1.80d0/





#--------------------------------------


# read xmol-type file
def readxmol(ifile,elem,xyz):
   """
   read xmol file
   """
   lines = ifile.readlines()
   nat = int(lines[0])
   title = lines[1]
   for l in lines[2:]:
       if l.split() ==[]:
         break
       type, x, y, z = l.split()
       xyz.append([float(x),float(y),float(z)])
       elem.append(type)
#       xyz.append(l)
   return nat

# write xmol-type file
def writexmol(name,nat,XYZ,frag):
   """
   write xmol file 
   """
   ofile = open( name, 'w')
   ofile.write(str(nat)+'\n')
   ofile.write(str(name)+'\n')
#   print >>ofile, str(nat)
#   print >>ofile, str(name)
   for i in frag[:]:
       ofile.write(str("% 5.5s % 4.12f % 4.12f % 4.12f \n" % (elem[i], float(XYZ[i,0]), float(XYZ[i,1]), float(XYZ[i,2]) )))
   ofile.close()
   return



def c_dist(di,dj): ##calculate distance between 2 lines of coords
        """
        cartesian distance between two vectors(coordinates). 
        """
        x=np.subtract(di,dj)
        dist=np.linalg.norm(x)
        return dist


def bond_mat(nat,elem,XYZ):
    """
    construct a bonding matrix (atom i, atom j). Bond is assumed when bond_length minus (cov_rad_i+cov_rad_j)/2
    is smaller then 0.5.
    """
    cov={'h': 0.6430, 'he': 0.6430,'li': 2.4570,'be': 1.9090,'b': 1.5870, 'c':1.4360,'n': 1.3090,\
       'o': 1.0960, 'f': 1.1200, 'ne': 0.9450, 'na': 2.9860,'mg': 2.6460,'al':2.4000,'si': 2.1920,\
       'p': 2.0600,'s': 1.8900,'cl': 1.7950,'ar': 1.7010,'k': 3.8360,'ca:' :3.2880,'sc':2.7210,\
       'ti': 2.4940, 'v': 2.3050, 'cr': 2.2300, 'mn': 2.2110,'fe': 2.2110,'co': 2.1920,'ni': 2.1730,\
       'cu': 2.2110,'zn': 2.3620, 'ga': 2.3810, 'ge': 2.3050, 'as': 2.2680,'se': 2.1920, 'br': 2.1540,\
       'kr': 2.1160,'rb': 4.0820, 'sr': 3.6090,'y': 3.0610,'zr': 2.7400,'nb': 2.5320,'mo': 2.4570,\
       'tc': 2.4000,'ru': 2.3620,'rh': 2.3620,'pd': 2.4190, 'ag': 2.5320, 'cd': 2.7970,'in': 2.7210,\
       'sn':  2.6650,'sb': 2.6460,'te': 2.5700,'i': 2.5130,'xe': 2.4760,'cs': 4.4410,'ba': 3.7420,'pb':2.740}
#       3.1940,3.1180,3.1180,3.0990,3.0800,3.0610,3.4960,
#       3.0420,3.0050,3.0050,2.9860,2.9670,2.9480,2.9480,
#       2.9480,2.7210,2.5320,2.4570,2.4190,2.3810,2.4000,
#       2.4570,2.5320,2.8160,2.7970,2.7780,2.7590,2.7590,
#       2.7400)

    
    bonds=[]
    for i in range(nat): 
        ei=str.lower(elem[i])
        for j in range(i+1,nat):
              ej=str.lower(elem[j])
              dist=c_dist(XYZ[i,:],XYZ[j,:])
              check=(float(cov[ei])+float(cov[ej]))*0.5
              if abs(dist-check) <= 0.5:
                   bonds.append((i,j))
    return bonds

def check_bond_lengths(bonds,XYZnew,XYZold,elem):
     status=0
     for i in bonds[:]:
        ai=i[0]
        aj=i[1]
        veci=XYZold[ai,:]
        vecj=XYZold[aj,:]
        distold=c_dist(veci,vecj)
        veca=XYZnew[ai,:]
        vecb=XYZnew[aj,:]
        distnew=c_dist(veca,vecb)
        if abs(distold-distnew) >= 0.01:
           print('ERROR in bond length: [atom1 atom2 delta_distance]', ai+1,'[',elem[ai],']',' - ',aj+1,'[',elem[aj],']',abs(distold-distnew))
           status=1
     return status



# --------------------------------------------------------------



def main():

  molname=args.molecule
  # read in coordinates
  f = open(molname, "r")
  nat = readxmol(f,elem,xyz)
  f.close()
  XYZ=np.array([xyz])
  XYZ.shape=(nat,3)

  XYZold=np.array(XYZ) # backup

  print( ' # atoms :',nat)
  #print ' requested operations :',' -> '.join(SYM[0])


  #set vars
  x1=args.cut[0]-1
  x2=args.cut[1]-1
  ax=(x1,x2)
  if x1 > -1:
     print( 'cutting bond:',x1+1,'[',elem[x1],']',' - ',x2+1,'[',elem[x2],']')


  # make bonding matrix
  bonds= bond_mat(nat,elem,XYZ)
  bondsOld=tuple(bonds) #backup

  # cut the bond 
  # requirement: x1<x2
  if x1 > -1:
    for b in bonds[:]:
      if ax == b: 
        bonds.remove(ax)

  tmp_list=[]
# care about isolated atoms
# a wee bit wonky...  ^(o_O)^
  for b in bonds[:]:
      tmp_list.append(b[0]) 
      tmp_list.append(b[1]) 
  test=[]
  lones=[]
  for i in range(0,nat):
       test=[]
       if i not in tmp_list:
         test.append(i)
         lones.append(test)  
  if len(lones) > 0:
     print( 'isolated atoms: ',lones)
  for i in range(0,len(lones)):
       print( 'atom'+str(i)+'.xyz', lones[i])
       writexmol('atom'+str(i)+'.xyz',1,XYZ,lones[i])
  
  print('')
  # process fragments
  # somehow we can end up with duplicates in the fragments, we remove them later with np.unique.
  mol=[0]
  frags=[]
  ifrag=np.zeros(maxfrag)
  found=1
  nr=0

  while bonds[:]:
    while found == 1:
      found=0
      for i in mol[:]:
        for j in bonds[:]:
          if i in j:
            if i == j[0]:
              mol.append(j[1])
            if i == j[1]:
              mol.append(j[0])
            bonds.remove(j)
            found=1
          atlist=np.unique(mol) #  removes duplicates!
         # print('frag:',nr,' : ', sorted(atlist))
    frags.append(mol)
    nr+=1
    if nr >=maxfrag+1:
      sys.exit("error: too many fragments found")
    if bonds[:]:
      mol=[bonds[0][0]]
      found=1
    else:
      break

       
  
  print('')
  print( 'writing fragment files')
  sum=0
  for f in range(0,nr):
     atlist=np.unique(frags[f]) #  removes duplicates!
     tnat=atlist.size
     sum=sum+tnat
     print( 'frag'+str(f)+'.xyz', "  -> #atoms ", tnat)
     writexmol('frag'+str(f)+'.xyz',tnat,XYZ,atlist)
  sum=sum+len(lones)
  if sum != nat:
    print( 'ERROR: some atoms are missing!')




if __name__ == '__main__':
  main()
