#include <iostream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Functions.h"






void GetResultsPO(int ne,Element netlist[MAX_ELEM],double solution[MAX_NODES+1]){
    int i;
    double i1,i2,ic,ie,ib,gbe,gbc,ibe,ibc,vb,vc,ve,isbe,isbc,vtbe,vtbc,
            ctbe,ctbc,c0be,c1be,c0bc,c1bc,alpha,alphaR,va;
    bool qtype;
    char c_qtype;

    i=1;
    while(i<=ne) {

        // Transistor

        if (netlist[i].getTransistorElement() != NULL) {
            string nameType = netlist[i].getName();
            Transistor* transistor;
            transistor = netlist[i].getTransistorElement();
            vc = solution[netlist[i].getA()];
            vb = solution[netlist[i].getB()];
            ve = solution[netlist[i].getC()];
            isbe = (*transistor).getIsBe();
            vtbe = (*transistor).getVtBe();
            isbc = (*transistor).getIsBc();
            vtbc = (*transistor).getVtBc();
            alpha = (*transistor).getAlpha();
            alphaR = (*transistor).getAlphaR();
            va = (*transistor).getVa();
            c0be = (*transistor).getC0Be();
            c1be = (*transistor).getC1Be();
            c0bc = (*transistor).getC0Bc();
            c1bc = (*transistor).getC1Bc();
            i1=isbe*(exp((vb-ve)/vtbe)-1);
            i2=isbc*(exp((vb-vc)/vtbc)-1);
            ic = (-i2+alpha*i1);
            //ic += (vc-ve)*ic/va;
            ie=(-i1+alphaR*i2);
            //(1-alphaR*alpha)
            ib = -(ic+ie);

            if (qtype==NPN)
                c_qtype = 'N';
            else
                c_qtype='P';
            ParametersNewtonRaphsonDiode(vb-ve, isbe, vtbe, &gbe, &ibe);
            ParametersNewtonRaphsonDiode(vb-vc, isbc, vtbc, &gbc, &ibc);
            DiodeResistanceandCapacitance(vb-ve,c0be,c1be,isbe, vtbe, &gbe,&ctbe);
            DiodeResistanceandCapacitance(vb-vc,c0bc,c1bc,isbc, vtbc, &gbc,&ctbc);
            cout << scientific<< nameType<<": tipo="<<c_qtype<<" Vbe="<<(vb-ve)<<" Vbc="<<(vb-vc)<<" Vce="<<(vc-ve)<<
                 " Ibe="<<i1<<" Ibc="<<i2<<" Ic="<<ic<<" Ib="<<ib<<" Cbe="<<ctbe<<" Cbc="<<ctbc<<endl;
        }
        i++;
    }

};




int SolveSystem(int nv,double Yn[MAX_NODES + 1][MAX_NODES + 2]) {
    int i, j, l, a;
    long double t, p;

    for (i = 1; i <= nv; i++) {
        t = 0.0;
        a = i;
        for (l = i; l <= nv; l++) {
            if (fabs(Yn[l][i]) > fabs(t)) {
                a = l;
                t = Yn[l][i];
            }
        }
        if (i != a) {
            for (l = 1; l <= nv + 1; l++) {
                p = Yn[i][l];
                Yn[i][l] = Yn[a][l];
                Yn[a][l] = (double)p;
            }
        }
        if (fabs(t) < TOLG) {
            cout << "Sistema singular\n";
            return 1;
        }
        for (j = nv + 1; j >= i; j--) {  /* Basta j>i em vez de j>0 */
            Yn[i][j] /= t;
            p = Yn[i][j];
            if (p != 0)  /* Evita operacoes com zero */
                for (l = 1; l <= nv; l++) {
                    if (l != i)
                        Yn[l][j] -= Yn[l][i] * p;
                }
        }
    }
    return 0;

}



/* Rotina que conta os nos e atribui numeros a eles */
int GetNodeNumber(char *name, int *nv, char lista[MAX_NODES + 1][MAX_NAME + 2] ) {
    int i, achou;

    i = 0;
    achou = 0;
    while (!achou && i <= *nv)
        if (!(achou = !strcmp(name, lista[i]))) i++;
    if (!achou) {
        if (*nv == MAX_NODES) {
            cout << "O programa so aceita ate" << *nv << "nos\n";
            exit(1);
        }
        (*nv)++;
        strcpy(lista[*nv], name);
        return *nv; /* novo no */
    } else {
        return i; /* no ja conhecido */
    };
}



void ParametersNewtonRaphsonDiode(double vn,double is,double vt, double *go,double*io){
    if (vn>V_MAX)
        vn = V_MAX;
    *go = (is/(vt))*exp(vn/(vt));
    *io = is*(exp(vn/(vt))-1)-vn*(*go);
}


void DiodeResistanceandCapacitance(double vn,double C0,double C1, double is, double vt, double *go, double *co){
    if (vn>100000)
        *go = (is/vt)*exp(V_MAX/vt);
    else *go = (is/vt)*exp(vn/vt);
    *co = 0;
    if (vn > 0) {
        *co += C1 * (exp(vn / vt) - 1);
    }
    if (vn > PHI/2) {
        *co += C0/sqrt(0.5);
    }
    else {
        *co += C0/sqrt(1- (vn/PHI));
    }

}


void ParametersNewtonRaphsonEarly(double vc, double vb, double ve,double gbe,double gbc, double ibe, double ibc,
                                  double alpha, double va, double *g1,double*g2,double*g3,double*io){
    *g1 = (alpha/va)*gbe*(vc - ve);
    //cout << "DEBUG G1 "<< *g1<<endl;
    *g2 = -gbc*(vc-ve)/va;
    //cout << "DEBUG G2 "<< *g2<<endl;
    *g3 = (alpha*(gbe*(vb-ve)+ibe)-(gbc*(vb-vc)+ibc))/va;
    //cout << "DEBUG G3 "<< *g3<<endl;
    *io = -((*g1) *(vb-ve)+(*g2)*(vb-vc));
    //cout << "DEBUG IO "<< *io<<endl;
}


int makeNetlistPO(Element netlist[MAX_ELEM], Element netlistPO[MAX_ELEM],int ne){
    int i,iPO,indexNR,new_ne;
    double value;
    char type;
    iPO =0;
    indexNR = 0;
    for (i=1;i<=ne;i++){
        iPO++;
        char nameType[MAX_NAME];
        strcpy(nameType,"");
        strcpy(nameType,netlist[i].getName());
        type = nameType[0];
        if (type=='C'){
            char new_name[MAX_NAME];
            strcpy(new_name,nameType);
            new_name[0]='R';
            strcat(new_name,"C");

            value = double(1)/TOLG;
            netlistPO[iPO].setName(new_name);
            netlistPO[iPO].setA(netlist[i].getA());
            netlistPO[iPO].setB(netlist[i].getB());
            netlistPO[iPO].setValue(value);
            netlistPO[iPO].setCapacitorElement(netlist[i].getCapacitorElement());
        }
        else if (type=='L'){

            char new_name[MAX_NAME];
            strcpy(new_name,nameType);
            new_name[0]='R';

            strcat(new_name,"L");
            value = TOLG;
            netlistPO[iPO].setName(new_name);
            netlistPO[iPO].setA(netlist[i].getA());
            netlistPO[iPO].setB(netlist[i].getB());
            netlistPO[iPO].setValue(value);
            netlistPO[iPO].setInductorElement(netlist[i].getInductorElement());
        }
        else if (type=='D'){
            double go, io;
            double is, vt;
            if (netlist[i].getDiodeElement()==NULL) {
                return INVALID_DIODE_ERROR;
            }

            Diode* diode;
            diode = netlist[i].getDiodeElement();
            is = (*diode).getIs();
            vt = (*diode).getVt();
            ParametersNewtonRaphsonDiode(INITIAL_VN_DIODE, is, vt, &go, &io);
            char new_Rname[MAX_NAME];
            char new_Iname[MAX_NAME];
            strcpy(new_Rname, nameType);
            strcpy(new_Iname, nameType);
            new_Rname[0] = 'R';
            new_Iname[0] = 'I';
            strcat(new_Rname,"D");
            strcat(new_Iname,"D");
            //value = TOLG;
            netlistPO[iPO].setName(new_Rname);
            netlistPO[iPO].setA(netlist[i].getA());
            netlistPO[iPO].setB(netlist[i].getB());
            netlistPO[iPO].setValue(double(1.0) / go);
            netlistPO[iPO].setDiodeElement(diode);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_Iname);
            netlistPO[iPO].setA(netlist[i].getA());
            netlistPO[iPO].setB(netlist[i].getB());
            netlistPO[iPO].setValue(io);
            netlistPO[iPO].setDiodeElement(diode);
            netlistPO[iPO].setNR(1);
            indexNR++;
        }
        else if (type=='Q'){

            double goBe, ioBe;
            double isBe, vtBe;
            double goBc, ioBc;
            double isBc, vtBc;
            double g1,g2,g3,io;
            double alpha,alphaR,va;
            bool qType;
            char new_RBename[MAX_NAME];
            char new_IBename[MAX_NAME];
            char new_RBcname[MAX_NAME];
            char new_IBcname[MAX_NAME];
            char new_GBename[MAX_NAME];
            char new_GBcname[MAX_NAME];
            char new_G1Cename[MAX_NAME];
            char new_G2Cename[MAX_NAME];
            char new_G3Cename[MAX_NAME];
            char new_ICename[MAX_NAME];
            char new_CBename[MAX_NAME];
            char new_CBcname[MAX_NAME];
            int sign;

            if (netlist[i].getTransistorElement()==NULL) {
                return INVALID_TRANSISTOR_ERROR;
            }

            Transistor* transistor;
            transistor = netlist[i].getTransistorElement();
            isBe = (*transistor).getIsBe();
            vtBe = (*transistor).getVtBe();
            isBc = (*transistor).getIsBc();
            vtBc = (*transistor).getVtBc();
            alpha = (*transistor).getAlpha();
            alphaR = (*transistor).getAlphaR();
            va = (*transistor).getVa();
            qType = (*transistor).getQType();
            sign = 1;
            if (qType==PNP)
                sign = -1;
            ParametersNewtonRaphsonDiode((INITIAL_VBQ-INITIAL_VEQ), isBe, vtBe, &goBe, &ioBe);
            ParametersNewtonRaphsonDiode((INITIAL_VBQ-INITIAL_VCQ), isBc, vtBc, &goBc, &ioBc);
            strcpy(new_RBename, nameType);
            strcpy(new_IBename, nameType);
            strcpy(new_RBcname, nameType);
            strcpy(new_IBcname, nameType);
            strcpy(new_GBename, nameType);
            strcpy(new_GBcname, nameType);
            strcpy(new_G1Cename, nameType);
            strcpy(new_G2Cename, nameType);
            strcpy(new_G3Cename, nameType);
            strcpy(new_CBename,nameType);
            strcpy(new_CBcname,nameType);
            strcpy(new_ICename, nameType);
            new_RBename[0] = 'R';
            new_IBename[0] = 'I';
            new_GBename[0] = 'G';
            new_CBename[0]='R';
            strcat(new_RBename,"QBE\0");
            strcat(new_IBename,"QBE\0");
            strcat(new_GBename,"QBE\0");
            strcat(new_CBename,"QBE\0");
            new_RBcname[0] = 'R';
            new_IBcname[0] = 'I';
            new_GBcname[0] = 'G';
            new_CBcname[0]='R';
            strcat(new_RBcname,"QBC\0");
            strcat(new_IBcname,"QBC\0");
            strcat(new_GBcname,"QBC\0");
            strcat(new_CBcname,"QBE\0");
            //Capacitores
            strcat(new_CBename,"C");
            strcat(new_CBcname,"C");
            netlistPO[iPO].setName(new_CBename);
            netlistPO[iPO].setA(netlist[i].getB());
            netlistPO[iPO].setB(netlist[i].getC());
            netlistPO[iPO].setValue(double(1)/TOLG);
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_CBcname);
            netlistPO[iPO].setA(netlist[i].getB());
            netlistPO[iPO].setB(netlist[i].getA());
            netlistPO[iPO].setValue(double(1)/TOLG);
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;

            netlistPO[iPO].setName(new_RBename);
            netlistPO[iPO].setA(netlist[i].getB()); // base
            netlistPO[iPO].setB(netlist[i].getC()); // emissor
            netlistPO[iPO].setValue(double(1.0) / goBe);
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_IBename);
            netlistPO[iPO].setA(netlist[i].getB()); // base
            netlistPO[iPO].setB(netlist[i].getC()); // emissor
            //netlistPO[iPO].setValue((ioBe));
            netlistPO[iPO].setValue((ioBe-alphaR*ioBc));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_RBcname);
            netlistPO[iPO].setA(netlist[i].getB()); // base
            netlistPO[iPO].setB(netlist[i].getA()); // coletor
            netlistPO[iPO].setValue(double(1.0) / goBc);
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_IBcname);
            netlistPO[iPO].setA(netlist[i].getB()); // base
            netlistPO[iPO].setB(netlist[i].getA()); // coletor
            //netlistPO[iPO].setValue((ioBc));
            netlistPO[iPO].setValue((ioBc-alpha*ioBe));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_GBename);
            netlistPO[iPO].setA(netlist[i].getC()); // emissor
            netlistPO[iPO].setB(netlist[i].getB()); // base
            netlistPO[iPO].setC(netlist[i].getB()); // base
            netlistPO[iPO].setD(netlist[i].getA()); // coletor
            netlistPO[iPO].setValue((alphaR*goBc));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_GBcname);
            netlistPO[iPO].setA(netlist[i].getA()); // coletor
            netlistPO[iPO].setB(netlist[i].getB()); // base
            netlistPO[iPO].setC(netlist[i].getB()); // base
            netlistPO[iPO].setD(netlist[i].getC()); // emissor
            netlistPO[iPO].setValue((alpha*goBe));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            // Efeito Early
            ParametersNewtonRaphsonEarly(INITIAL_VCQ, INITIAL_VBQ, INITIAL_VEQ,goBe,goBc,ioBe, ioBc,
                                         alpha, va, &g1,&g2,&g3,&io);
            new_G1Cename[0] = 'G';
            new_G2Cename[0] = 'G';
            new_G3Cename[0] = 'G';
            new_ICename[0] = 'I';
            strcat(new_G1Cename,"QCE1\0");
            strcat(new_G2Cename,"QCE2\0");
            strcat(new_G3Cename,"QCE3\0");
            strcat(new_ICename,"QCE\0");
            netlistPO[iPO].setName(new_G1Cename);
            netlistPO[iPO].setA(netlist[i].getA()); // coletor
            netlistPO[iPO].setB(netlist[i].getC()); // emissor
            netlistPO[iPO].setC(netlist[i].getB()); // base
            netlistPO[iPO].setD(netlist[i].getC()); // emissor
            netlistPO[iPO].setValue((g1));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_G2Cename);
            netlistPO[iPO].setA(netlist[i].getA()); // coletor
            netlistPO[iPO].setB(netlist[i].getC()); // emissor
            netlistPO[iPO].setC(netlist[i].getB()); // base
            netlistPO[iPO].setD(netlist[i].getA()); // coletor
            netlistPO[iPO].setValue((g2));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_G3Cename);
            netlistPO[iPO].setA(netlist[i].getA()); // coletor
            netlistPO[iPO].setB(netlist[i].getC()); // emissor
            netlistPO[iPO].setC(netlist[i].getA()); // coletor
            netlistPO[iPO].setD(netlist[i].getC()); // emissor
            netlistPO[iPO].setValue((g3));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);
            iPO++;
            netlistPO[iPO].setName(new_ICename);
            netlistPO[iPO].setA(netlist[i].getA()); // coletor
            netlistPO[iPO].setB(netlist[i].getC()); // emissor
            netlistPO[iPO].setValue((io));
            netlistPO[iPO].setTransistorElement(transistor);
            netlistPO[iPO].setNR(1);

            indexNR++;
        }
        else {
            netlistPO[iPO] = netlist[i];
        }
    }
    new_ne = iPO;

    return new_ne;
}


int  addCurrentAboveNode(int *nn, int *nv, int ne,Element netlist[MAX_ELEM],char lista[MAX_NODES + 1][MAX_NAME + 2]) {
    int i;
    char type;

    *nn = *nv;
    for (i = 1; i <= ne; i++) {
        string nameType = netlist[i].getName();
        type = nameType[0];
        if (type == 'V' || type == 'E' || type == 'F' || type == 'O') {
            (*nv)++;
            if (*nv > MAX_NODES) {
                return TOO_MANY_CURRENTS_ERROR;
                // "As correntes extra excederam o numero de variaveis permitido (" << MAX_NODES << ")\n";
            }
            strcpy(lista[*nv], "j"); /* Tem espaco para mais dois caracteres */
            strcat(lista[*nv], netlist[i].getName());
            netlist[i].setX(*nv);
        } else if (type == 'H') {
            *nv = *nv + 2;
            if (*nv > MAX_NODES) {
                return TOO_MANY_CURRENTS_ERROR;
            }
            strcpy(lista[*nv - 1], "jx");
            strcat(lista[*nv - 1], netlist[i].getName());
            netlist[i].setX(*nv - 1);
            strcpy(lista[*nv], "jy");
            strcat(lista[*nv], netlist[i].getName());
            netlist[i].setY(*nv);
        }
    }
    return 0;
}

void ListAll(int nv, int ne,Element netlist[MAX_ELEM],char lista[MAX_NODES + 1][MAX_NAME + 2]) {

    int i;
    char type;

    cout << "Variaveis internas: \n";
    for (i = 1; i <= nv; i++)
        cout << i << " -> " << lista[i] << endl;
    std::cin.get();
    cout << "Netlist interno final\n";
    for (i = 1; i <= ne; i++) {
        string nameType = netlist[i].getName();
        type = nameType[0];
        if (type == 'R' || type == 'I' || type == 'V') {

            cout << netlist[i].getName() << " " << netlist[i].getA() << " " << netlist[i].getB() << " "
                 << netlist[i].getValue() << endl;
        } else if (type == 'G' || type == 'E' || type == 'F' || type == 'H') {
            cout << netlist[i].getName() << " " << netlist[i].getA() << " " << netlist[i].getB() << " "
                 << netlist[i].getC() << " " << netlist[i].getD() << " " << netlist[i].getValue() << endl;
        } else if (type == 'O') {
            cout << netlist[i].getName() << " " << netlist[i].getA() << " " << netlist[i].getB() << " "
                 << netlist[i].getC() << " " << netlist[i].getD() << endl;
        }
        if (type == 'V' || type == 'E' || type == 'F' || type == 'O')
            cout << "Corrente jx: " << netlist[i].getX() << endl;
        else if (type == 'H')
            cout << "Correntes jx e jy: " << netlist[i].getX() << ", " << netlist[i].getY() << endl;
    }

}


int getLastNFromString(string original_string,int N, char* new_string){
    int i,size;
    size = original_string.length();
    i=0;
    if (size>=N) {
        for (i = size - N; i < size; i++) {
            new_string[i - size+ N] = original_string[i];
        }
        new_string[i-size+N]='\0';
    }
    else
        new_string[i]='\0';
    return 0;
}

void SetMNAPO(int nv,int ne,double Yn[MAX_NODES + 1][MAX_NODES + 2], Element netlist[MAX_ELEM],bool mode,bool NRmode,int indexNR){
    int i,j,k;
    char type;
    double g,vn,is,vt,go,io,g1,g2,g3,gbe,gbc,ibe,ibc,vb,vc,ve,isbe,isbc,vtbe,vtbc,alpha,alphaR,va;
    double previousSolution[MAX_NODES+1];
    char ID4[5],ID3[4]; // identificador de tamanho 3 ou 4 para identificar QBC/QBE/H1BC etc.. +1 \0


    /* Zera sistema */
    previousSolution[0]=0;
    for (i = 1; i <= nv; i++) {
        previousSolution[i]=Yn[i][nv+1];
        for (j = 0; j <= nv + 1; j++)
            Yn[i][j] = 0;
    }
    /* Monta estampas */
    for (i = 1; i <= ne; i++) {
        string nameType = netlist[i].getName();
        type = nameType[0];

        // Transistor

        if ((netlist[i].getTransistorElement()!=NULL)&&(NRmode==NR_ON)&&(mode==PO_MODE)){

            Transistor* transistor;
            transistor = netlist[i].getTransistorElement();
            /* Ativar newton Raphson */
            //cout << "DEBUG PS VC" << previousSolution[(*transistor).getA()]<<endl;
            //cout << "DEBUG PS VB" << previousSolution[(*transistor).getB()]<<endl;
            //cout << "DEBUG PS VE" << previousSolution[(*transistor).getC()]<<endl;
            vc = previousSolution[(*transistor).getA()];
            vb = previousSolution[(*transistor).getB()];
            ve = previousSolution[(*transistor).getC()];

            isbe = (*transistor).getIsBe();
            vtbe = (*transistor).getVtBe();
            isbc = (*transistor).getIsBc();
            vtbc = (*transistor).getVtBc();
            alpha = (*transistor).getAlpha();
            alphaR = (*transistor).getAlphaR();
            va = (*transistor).getVa();
            ParametersNewtonRaphsonDiode(vb-ve, isbe, vtbe, &gbe, &ibe);
            //cout << "DEBUG SetMNAPO GBE "<<gbe<<endl;
            //cout << "DEBUG SetMNAPO IBE "<<ibe<<endl;
            ParametersNewtonRaphsonDiode(vb-vc, isbc, vtbc, &gbc, &ibc);
            //cout << "DEBUG SetMNAPO GBC "<<gbc<<endl;
            //cout << "DEBUG SetMNAPO IBC "<<ibc<<endl;
            ParametersNewtonRaphsonEarly(vc, vb, ve,gbe,gbc,ibe, ibc,
                                         alpha, va, &g1,&g2,&g3,&io);

            if (type == 'R') {
                g=0;
                getLastNFromString(nameType,4,ID4);
                getLastNFromString(nameType,3,ID3);
                if (!strcmp(ID3,"QBE")){
                    g = gbe;
                }
                else if (!strcmp(ID3,"QBC")){
                    g=gbc;
                }
                else if (!strcmp(ID4,"QBEC")){
                    g=TOLG;
                }
                else if (!strcmp(ID4,"QBCC")){
                    g=TOLG;
                }

                Yn[netlist[i].getA()][netlist[i].getA()] += g;
                Yn[netlist[i].getB()][netlist[i].getB()] += g;
                Yn[netlist[i].getA()][netlist[i].getB()] -= g;
                Yn[netlist[i].getB()][netlist[i].getA()] -= g;
            } else if (type == 'I') {
                g=0;
                getLastNFromString(nameType,3,ID3);
                if (!strcmp(ID3,"QBE"))
                    g = ibe-alphaR*ibc;
                else if (!strcmp(ID3,"QBC")){
                    g=ibc-alpha*ibe;
                }
                else if (!strcmp(ID3,"QCE")){
                    g=io;
                }
                Yn[netlist[i].getA()][nv + 1] -= g;
                Yn[netlist[i].getB()][nv + 1] += g;
            } else if (type == 'G') {
                g=0;
                getLastNFromString(nameType,3,ID3);
                getLastNFromString(nameType,4,ID4);
                if (!strcmp(ID3,"QBE"))
                    g = alphaR*gbc;
                else if (!strcmp(ID3,"QBC"))
                    g=alpha*gbe;
                else if (!strcmp(ID4,"QCE1"))
                    g=g1;
                else if (!strcmp(ID4,"QCE2"))
                    g=g2;
                else if (!strcmp(ID4,"QCE3"))
                    g=g3;

                Yn[netlist[i].getA()][netlist[i].getC()] += g;
                Yn[netlist[i].getB()][netlist[i].getD()] += g;
                Yn[netlist[i].getA()][netlist[i].getD()] -= g;
                Yn[netlist[i].getB()][netlist[i].getC()] -= g;
            }
        }



            // Diode
        else if ((netlist[i].getDiodeElement()!=NULL)&&(NRmode==NR_ON)&&(mode==PO_MODE)){
            /* Ativar newton Raphson */
            vn = previousSolution[netlist[i].getA()] - previousSolution[netlist[i].getB()];

            Diode* diode;
            diode = netlist[i].getDiodeElement();
            is = (*diode).getIs();
            vt = (*diode).getVt();
            ParametersNewtonRaphsonDiode(vn, is, vt, &go, &io);
            if (type == 'R') {
                g = go;
                Yn[netlist[i].getA()][netlist[i].getA()] += g;
                Yn[netlist[i].getB()][netlist[i].getB()] += g;
                Yn[netlist[i].getA()][netlist[i].getB()] -= g;
                Yn[netlist[i].getB()][netlist[i].getA()] -= g;
            } else if (type == 'I') {
                g = io;
                Yn[netlist[i].getA()][nv + 1] -= g;
                Yn[netlist[i].getB()][nv + 1] += g;
            }
        }

        else if (type == 'R') {
            g = double(1) / netlist[i].getValue();
            Yn[netlist[i].getA()][netlist[i].getA()] += g;
            Yn[netlist[i].getB()][netlist[i].getB()] += g;
            Yn[netlist[i].getA()][netlist[i].getB()] -= g;
            Yn[netlist[i].getB()][netlist[i].getA()] -= g;
        } else if (type == 'G') {
            g = netlist[i].getValue();
            Yn[netlist[i].getA()][netlist[i].getC()] += g;
            Yn[netlist[i].getB()][netlist[i].getD()] += g;
            Yn[netlist[i].getA()][netlist[i].getD()] -= g;
            Yn[netlist[i].getB()][netlist[i].getC()] -= g;
        } else if (type == 'I') {
            g = netlist[i].getValue();
            Yn[netlist[i].getA()][nv + 1] -= g;
            Yn[netlist[i].getB()][nv + 1] += g;
        } else if (type == 'V') {
            Yn[netlist[i].getA()][netlist[i].getX()] += 1;
            Yn[netlist[i].getB()][netlist[i].getX()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getA()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getB()] += 1;
            Yn[netlist[i].getX()][nv + 1] -= netlist[i].getValue();
        } else if (type == 'E') {
            g = netlist[i].getValue();
            Yn[netlist[i].getA()][netlist[i].getX()] += 1;
            Yn[netlist[i].getB()][netlist[i].getX()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getA()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getB()] += 1;
            Yn[netlist[i].getX()][netlist[i].getC()] += g;
            Yn[netlist[i].getX()][netlist[i].getD()] -= g;
        } else if (type == 'F') {
            g = netlist[i].getValue();
            Yn[netlist[i].getA()][netlist[i].getX()] += g;
            Yn[netlist[i].getB()][netlist[i].getX()] -= g;
            Yn[netlist[i].getC()][netlist[i].getX()] += 1;
            Yn[netlist[i].getD()][netlist[i].getX()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getC()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getD()] += 1;
        } else if (type == 'H') {
            g = netlist[i].getValue();
            Yn[netlist[i].getA()][netlist[i].getY()] += 1;
            Yn[netlist[i].getB()][netlist[i].getY()] -= 1;
            Yn[netlist[i].getC()][netlist[i].getX()] += 1;
            Yn[netlist[i].getD()][netlist[i].getX()] -= 1;
            Yn[netlist[i].getY()][netlist[i].getA()] -= 1;
            Yn[netlist[i].getY()][netlist[i].getB()] += 1;
            Yn[netlist[i].getX()][netlist[i].getC()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getD()] += 1;
            Yn[netlist[i].getY()][netlist[i].getX()] += g;
        } else if (type == 'O') {
            Yn[netlist[i].getA()][netlist[i].getX()] += 1;
            Yn[netlist[i].getB()][netlist[i].getX()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getC()] += 1;
            Yn[netlist[i].getX()][netlist[i].getD()] -= 1;
        }
#ifdef DEBUG
        /* Opcional: Mostra o sistema apos a montagem da estampa */
        cout << "Sistema apos a estampa de " << netlist[i].getName() << " (n="<<indexNR<<")\n";
        for (k = 1; k <= nv; k++) {
            for (j = 1; j <= nv + 1; j++)
                if (Yn[k][j] != 0) //cout << Yn[k][j];
                    cout <<setw(11) << setprecision(2) << scientific << showpos << Yn[k][j]<<fixed;
                else cout << setw(11) << " ... ";
            cout << "\n";
        }
        //std::cin.get();
#endif
    }

}
