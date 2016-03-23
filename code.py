# -*- coding: cp1252 -*-
import random
from scipy import stats
from numpy import *
import matplotlib.pyplot as plt
import networkx as nx
import os
 
#-----------------------------------------------------------------------------------------------------------------------------
#            PARAMETRES
#-----------------------------------------------------------------------------------------------------------------------------
 
global TAUX1
global ELITISTE
global SEUIL
global GAMMA
 
#Taux d'arrette au départ
TAUX1=0.05
 
#Methode selection
ELITISTE=True
 
#nombre d'individu (=nombre de graphe)
N=150
 
#taux mutations
tauxMut=0.001
tauxCross=0.0002
 
#nombre de noeud par graphe
noeud=40
 
#Nb generation (SEUIL=False) ou seuil (SEUIL=True)
SEUIL=False
nbgen_ou_seuil=10
 
#Coeff petit monde
a1=1
 
#Coeff loi puissance
GAMMA=2.5  #entre 2 et 3
a2=1
a3=-0.2
 
#1-25: 0.3868  19.5564
 
 
#---------------------------------------------------------------------------------------------------------------------------
#         Classe INDIVIDU
#-----------------------------------------------------------------------------------------------------------------------------
class individu:
  #Constructeur d'un genome aleatoire ou constructeur de copie
  #nb : nombre de noeud, genome aleatoire de 0/1 avec TAUX1 de 1, fitness pour les deux prorietes
  def __init__(self,ind=None,nb=1):
    if ind is None :
      self.nb=nb
      self.genome=[]
      for i in range(sum(range(1,nb))):
        if random.random()<TAUX1:
          self.genome.append(1)
        else:
          self.genome.append(0)
    else:
      self.genome=ind.genome[:]
      self.nb=ind.nb


#A MODIFIER POUR LA CLIQUE
  #Calcul de la fitness de l'individu
  def fitness(self,a1,a2,a3):
    self.creation_graphe()

    f1=self.petit_monde()
    f2=self.loi_puissance()
    f3=self.loi_clustering()

    if f3==True:
      b=1
    else:
      b=0

    #on ne veut pas de graphe pas connecte
    #if nx.is_connected(self.graphe)==False:
     # sub=nx.connected_component_subgraphs(self.graphe)
      #f1*=len(sub)
      #f2*=len(sub)
    return a1*f1,a2*f2,a3*b 
 

    

  def petit_monde(self): 
   if nx.is_connected(self.graphe):
     l= nx.average_shortest_path_length(self.graphe)
     d=abs(l-log(log(self.nb)))
   else:
     sub=list(nx.connected_component_subgraphs(self.graphe))
     n=len(sub)
     l=0
     for g in sub:
       if g.number_of_nodes()==1:
         n=n-1
       else:
         l+=nx.average_shortest_path_length(g)
     l=l/1.0*n
     d=abs(l-log(log(self.nb)))*len(sub)
   return d
 

  def loi_puissance(self):
    distri = nx.degree_histogram(self.graphe)
    k_obs=nx.average_degree_connectivity(self.graphe)

    tab_deg=self.graphe.degree() #calcul les desgrés de tous les noeuds
    list_deg=[]
    type(tab_deg)
    for key,value in tab_deg.iteritems():
      temp = [key,value]
      list_deg.append(temp[1])    
 ########################""
    f=0
    for k in range(1,len(distri)):
      f+=abs(distri[k]-(self.nb* (k**(-GAMMA))))
 ########################
    pk=[x/self.nb for x in distri] 
    gradient, intercept, r_value, p_value, std_err = stats.linregress(distri,pk)
    #stats.kstest(k_theo,k_obs)

    return f



#CLIQUE

  def loi_clustering(self):

    #degre=nx.degree_centrality(self.graphe) 
    coef_k=nx.clustering(self.graphe,weight=None)
    moy_coef=nx.average_clustering(self.graphe)
    tab_deg=self.graphe.degree() #calcul les desgrés de tous les noeuds

    list_deg=[]
    
    list_deg_log=[]

    for key,value in tab_deg.iteritems():
      temp = [key,value]
      list_deg.append(temp[1])
      tmp=log(temp[1])
      list_deg_log.append(tmp)

    list_coef=[]
    list_coef_log=[]

    for key,value in coef_k.iteritems():
      temp = [key,value]
      list_coef.append(temp[1])
      tmp=log(temp[1])
      list_coef_log.append(tmp)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(list_deg_log,list_coef_log)

    distrib_th=-1*list_deg_log
    ks_val=stats.ks_2samp(list_coef_log,distrib_th)

    return ks_val>0.05


 #---------------------------------





  #Transformation d'un genome en graphe
  #genome : triangle supérieur de la matrice d'adjacence (sans la diagonale)
  def creation_graphe(self):
    self.graphe=nx.Graph()     #objet graphe
    # Création des noeuds
    self.graphe.add_nodes_from(range(self.nb))
    # Création des aretes
    l=self.nb-1
    i=0
    j=0
    k=0
    while l>0:
      for j in range(i+1,i+1+l):
        if self.genome[k]==1:
          self.graphe.add_edge(i,j)
        k+=1
      i+=1
      l-=1
 
 
#---------------------------------------------------------------------------------------------------------------------------------
#      Classe POPULATION
#---------------------------------------------------------------------------------------------------------------------------------
class pop:
  #population de N individu "type_ind", avec taux de mutation et taux de crossing
  def __init__(self,N,type_ind,txMut,txCross,nb):
    self.N=N                             #nombre d'individu (30,40 puis 100)
    self.txMut=txMut
    self.txCross=txCross
    self.type_ind=type_ind               #le type d'individu (ici on utilisera la classe "individu")
    self.pop=[type_ind(nb=nb) for i in range(N)]  #creation de la population : tableau des individus
    self.gen = 0                         #numero de generation
    self.tab_fitmoy=[[],[],[]]
    self.tab_fitmin=[[],[],[]]
 
  #calcul de la fitness de chaque individu et fitness moyenne/minimale de la population
  def calc_fitness(self,a1,a2,a3):
    self.f=[]        #tableau des fitness
    self.fitmoy=0     #fitness moyenne
    self.fitmin=[10000,0,0]  #fitness globale et detail de l'individu mmin
    self.petit_moy=0    #fitness petit monde moyenne
    self.loi_moy=0      #fitness loi puissance moyenne
    for i,x in enumerate(self.pop):  #i:indice  x:individu
      f1,f2,f3=x.fitness(a1,a2,a3)
      fi=f1+f2+f3
      self.fitmoy += fi
      self.petit_moy+=f1
      self.loi_moy+=f2
      if self.fitmin[0]>fi:
        self.fitmin=[fi,f1/a1,f2/a2]
      self.f.append([fi,i])  #fi:fitness, i:indice dans la population
    self.f.sort()            #tableau selon fitness croissante 
    self.fitmoy/=1.0*self.N
    self.petit_moy/=1.0*self.N
    self.loi_moy/=1.0*self.N
     
    self.tab_fitmoy[0].append(self.fitmoy)
    self.tab_fitmoy[1].append(self.petit_moy/a1)
    self.tab_fitmoy[2].append(self.loi_moy/a2)
     
    self.tab_fitmin[0].append(self.fitmin[0])
    self.tab_fitmin[1].append(self.fitmin[1])
    self.tab_fitmin[2].append(self.fitmin[2])
 
  #creation d'une population transitoire: on pioche N individus dans la population courante, en favorisant les plus basses fitness
  #l'index j permet de choisir un individu dans le tableau de fitness, en favorisant le début du tableau
  #self.f[j][1] donne l'indice de cet individu dans la population initiale
  def new_pop(self):
    self.npop=[]
    if not ELITISTE :
      for i in range(self.N):
        r= randint(0,(self.N+1)*(self.N)/2)
        ind = 0
        j=0
        while r>ind:
          ind += N-j
          j+=1
        self.npop.append(self.type_ind(self.pop[self.f[j-1][1]]))
    else:
      for i in range(10):
        self.npop.append(self.type_ind(self.pop[self.f[i][1]]))
      for i in range(self.N-10):
        r= random.randint(0,(self.N+1)*(self.N)/2)
        ind = 0
        j=0
        while r>ind:
          ind += N-j
          j+=1
        self.npop.append(self.type_ind(self.pop[self.f[j-1][1]]))
 
  #mutation de la population transitoire
  def mutation(self):
    for x in self.npop:
      g = x.genome
      for i in range(len(g)):
        if random.random()<self.txMut:
          g[i]=1-g[i]
    x.genome=g
 
  #crossing de la population transitoire
  #tous les genomes ont la même taille
  def cross(self):
    for x in self.npop:               #pour chaque individu
      if random.random()<self.txCross:
        g1 = x.genome                             #indivu 1 (parcourt de la population)
        g2= self.pop[random.randint(0,self.N-1)].genome  #individu 2 (choisit aléatoirement)
        z= random.randint(0,len(g1)-1)    #z : un endroit du genome de l'indivu 1
        if random.random()<0.5:
          g1[0:z]=g2[0:z]       #on copie genome 2 sur genome 1 (de 0 à z)
        else:
          g1[z:]=g2[z:]         #ou on copie genome 2 sur genome 1 (de z à la fin)
        x.genome=g1
 
  #Mise à jour : la population transitoire devient la population courante
  def update(self):
    self.pop=self.npop[:]
    self.gen+=1
 
  #Méthodes de sélection et de mutation
  def evolution(self):
    self.new_pop()
    self.mutation()
    self.cross()
 
  #Evolution d'une population sur X génération
  def loop(self,a1,a2,a3,nbgen_ou_seuil):
    if not SEUIL:
      while self.gen<=nbgen_ou_seuil:
        self.calc_fitness(a1,a2,a3)  #calcul tableau fitness, fitness min et moy 
        self.evolution()     #nouvelle population npop mutée
        self.update()        #pop=npop, generation +1
        print self.gen," fitmoy= ",round(self.fitmoy,4),"\t petit_moy= ",round(self.petit_moy/a1,4), "\t loi_moy= ",round(self.loi_moy/a2,4), "\t   individu min= ",round(self.fitmin[0],4),"\t",round(self.fitmin[1],4),"\t",round(self.fitmin[2],4)
    else :
      self.calc_fitness(a1,a2,a3)
      while self.fitmin>nbgen_ou_seuil:
        self.calc_fitness(a1,a2,a3)  #calcul tableau fitness, fitness min et moy 
        self.evolution()     #nouvelle population npop mutée
        self.update()        #pop=npop, generation +1
        print self.gen," fitmoy= ",round(self.fitmoy),"  fitmin= ",self.fitmin  #affichage de la génération, fitness moy et min
 
#-------------------------------------------------------------------------------------------------------------------------------
  #Sauvegarde les données de l'individu ayant la plus basse fitness
    # crée un fichier "nom".edge pour visualiser le graphe sous Cytoscape 
    # crée une image du graphe
    # crée un fichier texte contenant le génome de l'individu
  def save_Graph(self): # sauvegarde les données de l'individu ayant la plus basse fitness
    nbNode = self.pop[self.f[0][1]].nb
    nom="_N"+str(self.N)+"_Mut"+str(self.txMut)+"_Cross"+str(self.txCross)+"_noeud"+str(nbNode)+"_g"+str(self.gen)
    individu_min=self.pop[pop1.f[0][1]]
    individu_min.creation_graphe()    # on crée le graph
    G=individu_min.graphe # on récupère le graphe
     
    # création du fichier contenant le graghe
    tab = nx.adjacency_matrix(G) ## on crée la matrice d'adjacence
    fedge=open("data/graph"+ nom +".edges","w") # on ouvre le fichier de sauvegarde
    for i in range (0,nbNode) :
      for j in range (i,nbNode) :
        if tab[i,j] == 1:
          fedge.write(str(i+1)+"\t" + str(j+1)+ "\t" +" "+"\n") # enregistrement des edges (interactions entre noeud)
    fedge.close()
     
    # création de l'image (quand il y a trop de noeud devient ilisible)
    nx.draw(G)
    plt.savefig("data/image"+nom+".png") 
    plt.show()
    plt.close()
 
    # Enregiste le génome de l'individu ayant la plus basse fitness
    f=open("data/genome"+nom+".dat","w")
    s=""
    for i in individu_min.genome:
      s+=str(i)
    f.write(s)
    f.close()
    print "Genome : " + s
 
 
  def hist_Graph(self): # sauvegarde les données de l'individu ayant la plus basse fitness
    individu_min=self.pop[pop1.f[0][1]]
    distri = nx.degree_histogram(individu_min.graphe)
    liste = []
    print "distri"+ str(distri)
    for i in range(0,len(distri)) :
      for j in range(0,distri[i]) :
        liste.append(i)
    plt.hist(liste)
    plt.xlabel("degree")
    plt.ylabel("nombre de noeud")
    plt.show()
    plt.savefig("data/histogram.png")
    plt.close()
   
  def plot_Graph(self) :
   
    print str(len(range(1,nbgen_ou_seuil+1)))
    print str(len(self.tab_fitmin[0]))
   
    plt.plot(range(0,nbgen_ou_seuil+1),self.tab_fitmin[0],'b')
    plt.plot(range(0,nbgen_ou_seuil+1),self.tab_fitmoy[0],'r')
    plt.xlabel("nombre de generation")
    plt.ylabel("fitness")
    plt.show()
    plt.savefig("data/fitness.png")
    plt.close()
     
    plt.plot(range(0,nbgen_ou_seuil+1),self.tab_fitmin[1],'b')
    plt.plot(range(0,nbgen_ou_seuil+1),self.tab_fitmoy[1],'r')
    plt.xlabel("nombre de generation")
    plt.ylabel("coef f1")
    plt.show()
    plt.savefig("data/petitMonde.png")
    plt.close()
     
     
    plt.plot(range(0,nbgen_ou_seuil+1),self.tab_fitmin[2],'b')
    plt.plot(range(0,nbgen_ou_seuil+1),self.tab_fitmoy[2],'r')
    plt.xlabel("nombre de generation")
    plt.ylabel("coef f2")
    plt.show()
    plt.savefig("data/loipui.png")
    plt.close()
     
     
#MAIN
#---------------------------------------------------------------------------------------------------------------------
try:
  os.mkdir("data")
except OSError:
  pass
 
random.seed(11)
 
pop1=pop(N,individu,tauxMut,tauxCross,noeud)
pop1.loop(a1,a2,a3,nbgen_ou_seuil)
pop1.save_Graph()
pop1.hist_Graph()
pop1.plot_Graph()