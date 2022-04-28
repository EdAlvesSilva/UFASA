//
// Created by Carol Bez on 20-Jun-17.
//



#include <iostream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <complex>
#include "Functions.h"
#include "psFunctions.h"


int linspace(double begin,double end, double points,double frequencies[MAX_FREQUENCIES+1], int vb) {
    double step;

    step = (end - begin)/points;
    cout << "step : " << step << endl;
    int i=vb;
    if (i!=0) i++;
    for (;i<=points+vb;i++) {
        frequencies[i] = begin + (i-vb)*step;
        cout << setprecision(2)<< "indice i= " << i << "Recebeu : " << frequencies[i] << endl;
        cin.get();
    }
    return 0;
}

int makeFrequenciesList(double begin, double end, double points, double frequencies[MAX_FREQUENCIES + 1], double *number_of_frequencies, int type) {
    (*number_of_frequencies) = 0;
    if (type == LINEAR) {
        linspace(begin,end,points,frequencies);
        (*number_of_frequencies) = points;
        return 0;
    }

    if (type == LOG_DEC) {
        int i = 0;
        if (log10(end/begin) <= 1) {
            linspace(begin, end, points, frequencies);
            (*number_of_frequencies) = points;
            return 0;
        } else {
            double s;
            s = ceil(log10(begin));
            if (s != log10(begin)) {
                linspace(begin, pow(10, s), points, frequencies);
                i++;
                (*number_of_frequencies) += points;
            }
            cout << "Numero de frequencias : " << (*number_of_frequencies) << endl;
            cin.get();

            while (i <= (log10(end/begin) - 1)) {
                linspace(pow(10, (i + s - 1)), pow(10, i + s), points, frequencies, (int)(i * points));
                (*number_of_frequencies) += points;
                i++;
            }

            linspace(pow(10,(i+s-1)), end, points, frequencies, (int) (i * points));
            (*number_of_frequencies) += points;
            i++;
        }
        return 0;
    }

    if (type == LOG_OCT) {
        int i = 0;
        if((log(end/begin)/log(2))<= 1) {
            linspace(begin,end,points,frequencies);
            (*number_of_frequencies) = points;
            return 0;
        } else {
            double s;
            s = ceil(log(begin)/log(2));
            if (s != log(begin)/(log(2))) {
                linspace(begin, pow(2, s), points, frequencies);
                i++;
                (*number_of_frequencies) += points;
            }

            while (i <= ((log(end / begin)/log(2)) - 1)) {
                linspace(pow(2, i + s-1), pow(2, i + s ), points, frequencies, (int) (i * points));
                (*number_of_frequencies) += points;
                i++;
            }

            linspace(pow(2,i+s-1), end, points, frequencies, (int) (i * points));
            (*number_of_frequencies) += points;
            i++;
            return 0;
        }
    }
    else {
        return NOT_ANAL_TYPE_ERROR;
    }
}

int makeNetlistPS(Element netlist[MAX_ELEM], Element netlistPS[MAX_ELEM], int ne, double solution[MAX_NODES + 1]){
    int i,iPS,new_ne;
    double c0,c1, vn;
    char type;
    iPS =0;


    for (i=1;i<=ne;i++){
        iPS++;
        char nameType[MAX_NAME];
        strcpy(nameType,"");
        strcpy(nameType,netlist[i].getName());
        type = nameType[0];
        if (type=='D'){
            double go, co;
            double is, vt;
            if (netlist[i].getDiodeElement()==NULL) {
                return INVALID_DIODE_ERROR;
            }

            Diode* diode;
            diode = netlist[i].getDiodeElement();
            c0 = (*diode).getC0();
            c1 = (*diode).getC1();
            vt = (*diode).getVt();
            is = (*diode).getIs();
            vn = solution[netlist[i].getA()] - solution[netlist[i].getB()];
            if (vn>V_MAX) vn = V_MAX;
            DiodeResistanceandCapacitance(vn,c0,c1, is, vt, &go, &co);
            char new_Rname[MAX_NAME];
            char new_Cname[MAX_NAME];
            strcpy(new_Rname, nameType);
            strcpy(new_Cname, nameType);
            new_Rname[0] = 'R';
            new_Cname[0] = 'C';
            netlistPS[iPS].setName(new_Rname);
            netlistPS[iPS].setA(netlist[i].getA());
            netlistPS[iPS].setB(netlist[i].getB());
            netlistPS[iPS].setValue(1 / go);
            netlistPS[iPS].setDiodeElement(diode);
            netlistPS[iPS].setNR(0);
            iPS++;
            netlistPS[iPS].setName(new_Cname);
            netlistPS[iPS].setA(netlist[i].getA());
            netlistPS[iPS].setB(netlist[i].getB());
            netlistPS[iPS].setValue(co);
            netlistPS[iPS].setDiodeElement(diode);
        }
        if (type == 'Q') {
            Transistor* transistor;
            transistor = netlist[i].getTransistorElement();
            double vbe, vce,vbc, cd, cr,c0d,c1d,c0r,c1r,isd,isr, gd, gr, va, alpha_d, alpha_r,vtd,vtr;
            bool qType;
            int sign = 1;
            vbe = solution[netlist[i].getB()] - solution[netlist[i].getC()];
            vce = solution[netlist[i].getA()] - solution[netlist[i].getC()];
            vbc = solution[netlist[i].getB()] - solution[netlist[i].getA()];
            qType = (*transistor).getQType();
            isd = (*transistor).getIsBe();
            isr = (*transistor).getIsBc();
            c0d = (*transistor).getC0Be();
            c1d = (*transistor).getC1Be();
            c0r = (*transistor).getC0Bc();
            c1r = (*transistor).getC1Bc();
            vtd = (*transistor).getVtBe();
            vtr = (*transistor).getVtBc();
            alpha_d = (*transistor).getAlpha();
            alpha_r = (*transistor).getAlphaR();
            va = (*transistor).getVa();
            if (qType == PNP) sign = -1;
            DiodeResistanceandCapacitance(vbe,c0d,c1d,isd,vtd,&gd,&cd);
            DiodeResistanceandCapacitance(vbc,c0r,c1r,isr,vtr,&gr,&cr);

            char new_RCename[MAX_NAME];
            char new_RBename[MAX_NAME];
            char new_RBcname[MAX_NAME];
            char new_GBename[MAX_NAME];
            char new_GBcname[MAX_NAME];
            char new_CBename[MAX_NAME];
            char new_CBcname[MAX_NAME];

            strcpy(new_RCename, nameType);
            strcpy(new_RBename, nameType);
            strcpy(new_RBcname, nameType);
            strcpy(new_GBename, nameType);
            strcpy(new_GBcname, nameType);
            strcpy(new_CBename,nameType);
            strcpy(new_CBcname,nameType);
            new_RCename[0] = 'R';
            new_RBename[0] = 'R';
            new_GBename[0] = 'G';
            new_CBename[0]='C';
            strcat(new_RCename,"QCE\0");
            strcat(new_RBename,"QBE\0");
            strcat(new_GBename,"QBE\0");
            strcat(new_CBename,"QBE\0");
            new_RBcname[0] = 'R';
            new_GBcname[0] = 'G';
            new_CBcname[0]='C';
            strcat(new_RBcname,"QBC\0");
            strcat(new_GBcname,"QBC\0");
            strcat(new_CBcname,"QBC\0");
            //Capacitores
            strcat(new_CBename,"C");
            strcat(new_CBcname,"C");

            netlistPS[iPS].setName(new_CBename);
            netlistPS[iPS].setA(netlist[i].getB());
            netlistPS[iPS].setB(netlist[i].getC());
            netlistPS[iPS].setValue(cd);

            iPS++;

            netlistPS[iPS].setName(new_CBcname);
            netlistPS[iPS].setA(netlist[i].getB());
            netlistPS[iPS].setB(netlist[i].getA());
            netlistPS[iPS].setValue(cr);
            iPS++;

            netlistPS[iPS].setName(new_RBename);
            netlistPS[iPS].setA(netlist[i].getB()); // base
            netlistPS[iPS].setB(netlist[i].getC()); // emissor
            netlistPS[iPS].setValue(double(1.0)/gd);
            iPS++;

            netlistPS[iPS].setName(new_RBcname);
            netlistPS[iPS].setA(netlist[i].getB()); // base
            netlistPS[iPS].setB(netlist[i].getA()); // coletor
            netlistPS[iPS].setValue(double(1.0) / gr);
            iPS++;

            netlistPS[iPS].setName(new_GBename);
            netlistPS[iPS].setA(netlist[i].getC()); // emissor
            netlistPS[iPS].setB(netlist[i].getB()); // base
            netlistPS[iPS].setC(netlist[i].getB()); // base
            netlistPS[iPS].setD(netlist[i].getA()); // coletor
            netlistPS[iPS].setValue((alpha_r*gr));
            iPS++;

            netlistPS[iPS].setName(new_GBcname);
            netlistPS[iPS].setA(netlist[i].getA()); // coletor
            netlistPS[iPS].setB(netlist[i].getB()); // base
            netlistPS[iPS].setC(netlist[i].getB()); // base
            netlistPS[iPS].setD(netlist[i].getC()); // emissor
            netlistPS[iPS].setValue((alpha_d*gd));
            iPS++;

            //efeito early

            double i1, i2;
            if (qType == NPN) {

                i1=isd*(exp((vbe)/vtd)-1);
                i2=isr*(exp((vbc)/vtr)-1);


            }

            netlistPS[iPS].setName(new_RCename);
            netlistPS[iPS].setA(netlist[i].getA()); // coletor
            netlistPS[iPS].setB(netlist[i].getC()); // emissor
            netlistPS[iPS].setValue(va/(-i2+alpha_d*i1));

        }
        else {
            netlistPS[iPS] = netlist[i];
        }
    }
    new_ne = iPS;

    for (i=1;i<=new_ne;i++) {
        cout << netlistPS[i].getValue() << endl;
    }
    return new_ne;
}

void SetMNAPS(int nv,int ne, std::complex<double> Yn[MAX_NODES + 1][MAX_NODES + 2], Element netlist[MAX_ELEM],double frequency){
    int i,j;
    char type;
    double g,l;
    std::complex<double> temp;
    /* Zera sistema */
    for (i = 0; i <= nv+1; i++) {
        for (j = 0; j <= nv + 1; j++)
            Yn[i][j] = 0;
    }
    /* Monta estampas */
    for (i = 1; i <= ne; i++) {
        string nameType = netlist[i].getName();
        type = nameType[0];
        if (type == 'C'){
            g = (netlist[i].getValue() * 2 * PI) * frequency;
            Yn[netlist[i].getA()][netlist[i].getA()] += std::complex<double>(0,g);
            Yn[netlist[i].getB()][netlist[i].getB()] += std::complex<double>(0,g);
            Yn[netlist[i].getA()][netlist[i].getB()] -= std::complex<double>(0,g);
            Yn[netlist[i].getB()][netlist[i].getA()] -= std::complex<double>(0,g);
        }
        if (type == 'L'){
            g = -double(1)/((netlist[i].getValue())*PI*2*frequency);
            Yn[netlist[i].getA()][netlist[i].getA()] += std::complex<double>(0,g);
            Yn[netlist[i].getB()][netlist[i].getB()] += std::complex<double>(0,g);
            Yn[netlist[i].getA()][netlist[i].getB()] -= std::complex<double>(0,g);
            Yn[netlist[i].getB()][netlist[i].getA()] -= std::complex<double>(0,g);
        }
        if (type == 'R') {
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
            Source* source;
            source = netlist[i].getSourceElement();
            g = (*source).getModule();
            l = (*source).getPhase();
            temp = std::complex<double>(g*cos(l),g*sin(l));
            Yn[netlist[i].getA()][nv + 1] -= temp;
            Yn[netlist[i].getB()][nv + 1] += temp;
        } else if (type == 'V') {
            Source* source;
            source = netlist[i].getSourceElement();
            g = (*source).getModule();
            l = (*source).getPhase();
            temp = std::complex<double>(g*cos(l),g*sin(l));
            Yn[netlist[i].getA()][netlist[i].getX()] += 1;
            Yn[netlist[i].getB()][netlist[i].getX()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getA()] -= 1;
            Yn[netlist[i].getX()][netlist[i].getB()] += 1;
            Yn[netlist[i].getX()][nv + 1] -= temp;
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
        cout << "Sistema apos a estampa de " << netlist[i].getName() << endl;
        int k;
        for (k = 1; k <= nv; k++) {
            for (j = 1; j <= nv + 1; j++)
                if (abs(Yn[k][j]) != 0) cout << scientific << setprecision(3) <<"[" <<Yn[k][j].real()<< "  j" << Yn[k][j].imag() << "]";
                else cout << " ... ";
            cout << "\n";
        }
        //std::cin.get();
#endif
    }

}


