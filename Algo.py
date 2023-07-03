import numpy as np
dimension=5;



def creation_a(array_2d) :
    array_2d[0, 0] = array_2d[1, 0] = array_2d[dimension, 0] = 1;
    array_2d[dimension * dimension - 1, dimension * dimension - 1] = 1;
    array_2d[dimension * dimension - 2, dimension * dimension - 1] = 1;
    array_2d[dimension * (dimension - 1) - 1, (dimension * dimension - 1)] = 1;

    for k in range(1,dimension*dimension-1) : #de 1 a 23 si 5x5
        array_2d[k,k] = 1;
        if (k+1)//dimension == k//dimension :
            array_2d[k+1, k] = 1
        else :
            array_2d[k + dimension, k ] = 1;

        if (k-1)//dimension == k//dimension :
            array_2d[k-1, k] = 1;
        else :
            array_2d[k - dimension, k] = 1;

        if k+dimension <= dimension*dimension-1 :
            array_2d[k+dimension, k] = 1;

        if k-dimension >= 0 :
            array_2d[k-dimension, k] = 1;

#===========================================================
#===========================================================
#===========================================================
#===========================================================
#===========================================================
#===========================================================
#===========================================================

def MoinsZ_2Z(XA,XB) :
    if XA == XB :
        return 0;
    else :
        return 1;

def LaMoinsLambdaLb(A,La,Lb,Lambda):
    """Applique l’opération élémentaire La - lamda * Lb à la matrice A"""
    n = len(A[0]) # nbre de colonnes de A
    for k in range(n):
       A[La][k] = MoinsZ_2Z(A[La][k],Lambda*A[Lb][k]);
    return A;

def TriangulationSup(A,B) :
    """transforme A en matrice triangulaire supérieure (A doit avoir ça diagonale remplie)"""
    coef = 0;
    for i in range(0,len(A)): # pour chaque colonnes
        for j in range(i+1,len(A)) : # Pour les [i,25] dernières lignes | A[2,1] = colone 1 ligne 2
            if A[i,i] == 0 :
                k=i+1
                leave = False;
                while A[k,i] == 0 and leave == False : # trouve un ligne ou le coef de la colonne i n'est pas nul
                    if k+1 >= len(A) : # si on a pas trouver de coef on sort
                        leave = True;
                    else :
                        k = k+1;

                if A[k,i] != 0 : # si on a trouver coef alors on interverti les lignes
                    IntervertirLigne(A,i,k);
                    IntervertirLigne(B,i,k);

            # ========================================================================================
            if A[i,i] != 0 :
                coef = A[j,i] / A[i,i];
                LaMoinsLambdaLb(A,j,i,coef);
                LaMoinsLambdaLb(B, j, i, coef);
                #print("Ligne " + str(j) + " = Ligne " + str(j) + " - " + str(coef) + " * Ligne" + str(i));
                #print(A);
                #print("==========Sur Diag==========");
                #print(B);
                # ========================================================================================



def IntervertirLigne(A,La,Lb) :
    for i in range(0,len(A)) :
        Aux = A[La,i];
        A[La,i] = A[Lb,i];
        A[Lb,i] = Aux;
    #print("Ligne " + str(La) + " et Ligne " + str(Lb) + " On été interchanger ");
    #print(A);



def Diagonalisation(A,B) :
    """transforme A en matrice diagonale (A doit avoir ça diagonale remplie)"""
    coef =0;
    for i in range(len(A)-1,-1,-1): # pour chaque colonnes en partant de la fin
        for j in range(i-1,-1,-1) : # Pour les i premières lignes en partant de la dernière
            if A[i, i] != 0:
                coef = A[j,i] / A[i,i];
                LaMoinsLambdaLb(A,j,i,coef);
                LaMoinsLambdaLb(B,j,i,coef);
                #print("Ligne " + str(j) + " = Ligne " + str(j) + " - " + str(coef) + " * Ligne" + str(i));
                #print(A);
                #print("==========Sur Diag==========");
                #print(B);

import sys
def sortie(A):
    with open('Output.txt', 'w') as f:
        for i in range(len(A)):
            output = "";
            for j in range(len(A)):
                output = output + str(round(A[i,j],2)) + " , ";
            print(output, file=f);

        for k in range(len(A)):
            if A[k,k] == 0 :
                output = "Vecteur u" + str(k) + " = ";
                for h in range(len(A)) :
                    output = output + str(round(A[h,k],2)) + " , ";
                print(output, file=f);


def NormalisationDeLaDiag(A,B) : # B = matrice I à l'origine
    for i in range(0,len(A)): # pour chaque lignes
        aux = A[i, i];
        for j in range(0, len(A)): # pour chaque colonne
            if aux != 0: # si le coef dans la diag non nul
                B[i, j] = B[i, j] / aux;
                A[i, j] = A[i, j] / aux;


def Solution(R,B) :
    Out = np.zeros((len(R), 1), order='F');
    for i in range(len(R)) :
        for j in range(len(R)) :
            Out[i] = MoinsZ_2Z(Out[i],R[i,j] * B[j]);
    return Out;

def Optimale(A,Out) :
    B = Out; #on initialise Out comme solution optimale*
    nomOpti = "B";
    nom = ""
    print("la solution B +" + nom + "prend " + str(TotalDeCoup(B)) + "coups");
    for k in range(len(A)):
        if A[k,k] == 0 : #S'il s'agit d'une action nulle
            output = np.zeros((len(A), 1), order='F');
            for h in range(len(A)) :
                output[h] = A[h,k];
            output[k] = 1;
            nom = "U" + str(k);

            #calcul de C qui est la solution suivante
            C = np.zeros((len(A), 1), order='F');
            for l in range(len(A)) :
                C[l] = MoinsZ_2Z(B[l],output[l]); #C = B[l] + vecteur u... on obtient une new solution

            print("la solution B +" + nom + "prend " + str(TotalDeCoup(C)) + "coups");
            if TotalDeCoup(B) > TotalDeCoup(C) : #si la solution optimale est plus longue que la solution C
                B = C;
                nomOpti = "B " + nom;
    print("la solution optimale est donc " + nomOpti);
    s="";
    for i in range(len(B)) :
        s = s + str(int(B[i])) + ","
    print(nomOpti + " = " + s);

def TotalDeCoup(A) :
    output = 0;
    for i in range(len(A)) :
        if A[i] == 1 :
            output = output + 1;
    return output;

# ===========================================================
# ===========================================================
# ===========================================================
# ===========================================================
# ===========================================================
# ===========================================================
# ===========================================================

tableau = np.zeros((dimension*dimension,dimension*dimension),order='F'); # 0.0 a 24.24
Diag = np.zeros((dimension*dimension,dimension*dimension),order='F');

M = np.zeros((dimension*dimension,1),order='F');
M = [1,0,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,1,0,1,1,1,1,0,1];

# tableau2 = np.zeros((dimension*dimension,dimension*dimension),order='F');
tableau2 = np.matrix([[1.0,0.0,2.0,4.0],[0.0,1.0,3.0,2.0],[2.0,-1.0,0.0,3.0],[4.0,-1.0,3.0,9.0]]);

np.fill_diagonal(Diag,1);
creation_a(tableau);

cible=tableau;


#print(cible);
#print("=====================");
#print(Diag);

TriangulationSup(cible,Diag);
print("===========================================================");
print("=============== Debut de la Diagonalisation ===============");
print("===========================================================");

Diagonalisation(cible,Diag);

NormalisationDeLaDiag(cible,Diag);

#print("===========================================================");
#print("====================== Resultat Final =====================");
#print("======================== Matrice  A =======================");

#print(cible);
#print("===================== Matrice  Inverse ====================");
Optimale(cible,Solution(Diag,M));

sortie(cible);






#print(tableau2);