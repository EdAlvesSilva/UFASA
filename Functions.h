//
// Created by Carol Bez on 16-Jun-17.
//

#ifndef UNTITLED_FUNCTIONS_H
#define UNTITLED_FUNCTIONS_H




#define MAX_LINE 90
#define MAX_NAME 11
#define MAX_ELEM 50
#define MAX_NODES 50
#define TOLG 1e-9
#define DEBUG
#define PO_MODE 0
#define PS_MODE 1
#define NR_ON 1
#define NR_OFF 0
#define NPN 0
#define PNP 1
#define NO_SUCH_FILE_ERROR 1
#define INVALID_DIODE_ERROR 2
#define TOO_MANY_CURRENTS_ERROR 3
#define TRANSISTOR_INVALID_TYPE_ERROR 4
#define INVALID_TRANSISTOR_ERROR 5
#define NOT_ANAL_TYPE_ERROR   6



// Default Values
#define START_POINT 0
#define END_POINT 1000
#define NUMBER_OF_POINTS 100
#define LINEAR 0
#define LOG_DEC 1
#define LOG_OCT 2
#define PHI 0.6
#define N_VALUE 0.5
#define CO_DEF 5e-12
#define CI_DEF 100e-18
#define IS_DEF 1e-9
#define VT_DEF 43.42944e-3
#define ALPHA_DEF 0.99
#define ALPHAR_DEF 0.5
#define VA_DEF 100


#define INITIAL_VN_DIODE 0.1
#define INITIAL_VCQ 0.1
#define INITIAL_VBQ 0.1
#define INITIAL_VEQ 0.1
#define V_MAX 1.363



using namespace std;





class Source;
class Diode;
class Transistor;
class Capacitor;
class Inductor;
class Coupling;




class Element {
public:

    Element() {
        diodeElement = NULL;
        capacitorElement = NULL;
        inductorElement = NULL;
        sourceElement = NULL;
        couplingElement = NULL;
        transistorElement = NULL;
    }


    void setName(char val_name[MAX_NAME]) {
        strcpy(name, val_name);
    }

    char *getName() {
        return name;
    }

    void setA(int a_val) {
        a = a_val;
    }

    int getA() {
        return a;
    }

    void setB(int b_val) {
        b = b_val;
    }

    int getB() {
        return b;
    }

    void setC(int c_val) {
        c = c_val;
    }

    int getC() {
        return c;
    }

    void setD(int d_val) {
        d = d_val;
    }

    int getD() {
        return d;
    }

    void setX(int x_val) {
        x = x_val;
    }

    int getX() {
        return x;
    }

    void setY(int y_val) {
        y = y_val;
    }

    int getY() {
        return y;
    }

    void setValue(double val_value) {
        value = val_value;
    }

    double getValue() {
        return value;
    }


    void setNR(bool val_NR) {
        NR=val_NR;
    }

    char getNR() {
        return NR;
    }

    void setDiodeElement(Diode* val_diodeElement){
        diodeElement = val_diodeElement;
    }
    Diode* getDiodeElement(){
        return diodeElement;
    }

    void setSourceElement(Source* val_sourceElement){
        sourceElement = val_sourceElement;
    }
    Source* getSourceElement(){
        return sourceElement;
    }

    void setCapacitorElement(Capacitor* val_capacitorElement){
        capacitorElement = val_capacitorElement;
    }
    Capacitor* getCapacitorElement(){
        return capacitorElement;
    }

    void setInductorElement(Inductor* val_inductorElement){
        inductorElement = val_inductorElement;
    }
    Inductor* getInductorElement(){
        return inductorElement;
    }

    void setCouplingElement(Coupling* val_couplingElement){
        couplingElement = val_couplingElement;
    }
    Coupling* getCouplingElement(){
        return couplingElement;
    }
    void setTransistorElement(Transistor* val_transistorElement){
        transistorElement = val_transistorElement;
    }
    Transistor* getTransistorElement(){
        return transistorElement;
    }

protected:

    Diode* diodeElement;
    Capacitor* capacitorElement;
    Inductor* inductorElement;
    Coupling* couplingElement;
    Source* sourceElement;
    Transistor* transistorElement;
    char name[MAX_NAME];
    double value;
    int a, b, c, d, x, y;
    bool NR; /* if NR==1, value may change in Newton Raphson */


};









class Source {
public:
    Source(){};
    void setModule(double val_mod){
        module = val_mod;
    }
    double getModule (){
        return module;
    }
    void setPhase (double val_mod){
        phase = val_mod;
    }
    double getPhase (){
        return phase;
    }

    void setValue(double val_value) {
        value = val_value;
    }

    double getValue() {
        return value;
    }
protected:
    double module, phase,value;
};


class Diode : public Element {
public:

    Diode(){
    };
    void setIs(double val_is) {
        is = val_is;
    }

    double getIs() {
        return is;
    }
    void setVt(double val_vt) {
        vt = val_vt;
    }

    double getVt() {
        return vt;
    }

    void setC0(double val_c0) {
        c0 = val_c0;
    }

    double getC0() {
        return c0;
    }

    void setC1(double val_c1) {
        c1 = val_c1;
    }

    double getC1() {
        return c1;
    }
protected:
    double is,vt,c0,c1;

};


class Transistor : public Element{
public:

    Transistor(){
        Element();
    };
    void setIsBe(double val_isBe) {
        isBe = val_isBe;
    }

    double getIsBe() {
        return isBe;
    }
    void setVtBe(double val_vtBe) {
        vtBe = val_vtBe;
    }

    double getVtBe() {
        return vtBe;
    }

    void setC0Be(double val_c0Be) {
        c0Be = val_c0Be;
    }

    double getC0Be() {
        return c0Be;
    }

    void setC1Be(double val_c1Be) {
        c1Be = val_c1Be;
    }

    double getC1Be() {
        return c1Be;
    }
    void setIsBc(double val_isBc) {
        isBc = val_isBc;
    }

    double getIsBc() {
        return isBc;
    }
    void setAlpha(double val_alpha) {
        alpha = val_alpha;
    }

    double getAlpha() {
        return alpha;
    }
    void setVa(double val_va) {
        va = val_va;
    }

    double getVa() {
        return va;
    }


    void setAlphaR(double val_alphaR) {
        alphaR = val_alphaR;
    }

    double getAlphaR() {
        return alphaR;
    }

    void setVtBc (double val_vtBc) {
        vtBc = val_vtBc;
    }

    double getVtBc() {
        return vtBc;
    }

    void setC0Bc(double val_c0Bc) {
        c0Bc = val_c0Bc;
    }

    double getC0Bc() {
        return c0Bc;
    }

    void setC1Bc(double val_c1Bc) {
        c1Bc = val_c1Bc;
    }

    double getC1Bc() {
        return c1Bc;
    }

    void setQType(bool val_qType) {
        qType = val_qType;
    }

    bool getQType() {
        return qType;
    }




protected:
    double isBe,vtBe,c0Be,c1Be,isBc,vtBc,c0Bc,c1Bc,va,alpha,alphaR;
    bool qType;

};



class Capacitor {
public:
    Capacitor(){};
protected:
    bool mode;
/* mode = NR_OFF para analise do ponto de operacao ou 1 para resposta em frequencia */

};

class Inductor {
public:
    Inductor(){};
    void setCouple (bool coupling){
        couple = coupling;
    }
    bool checkCouple() {
        return couple;
    }
protected:
    bool mode;
/* mode = NR_OFF para analise do ponto de operacao ou 1 para resposta em frequencia */
    bool couple;
};


class Coupling {
public:
    Coupling(){};
    void setNameLA(char val_name[MAX_NAME]) {
        strcpy(nameLA, val_name);
    }

    char *getNameLA() {
        return nameLA;
    }
    void setNameLB(char val_name[MAX_NAME]) {
        strcpy(nameLB, val_name);
    }

    char *getNameLB() {
        return nameLB;
    }
    void setValue(double val_value) {
        value = val_value;
    }

    double getValue() {
        return value;
    }

protected:
    char nameLA[MAX_NAME];
    char nameLB[MAX_NAME];
    double value;
};






template <typename type>
int SolveSystem(int nv,type Yn[MAX_NODES + 1][MAX_NODES + 2]) {
    int i, j, l, a;
    type t, p;

    for (i = 1; i <= nv; i++) {
        t = 0.0;
        a = i;
        for (l = i; l <= nv; l++) {
            if (abs(Yn[l][i]) > abs(t)) {
                a = l;
                t = Yn[l][i];
            }
        }
        if (i != a) {
            for (l = 1; l <= nv + 1; l++) {
                p = Yn[i][l];
                Yn[i][l] = Yn[a][l];
                Yn[a][l] = p;
            }
        }
        if (abs(t) < TOLG) {
            cout << "Sistema singular\n";
            return 1;
        }
        for (j = nv + 1; j >= i; j--) {  /* Basta j>i em vez de j>0 */
            Yn[i][j] /= t;
            p = Yn[i][j];
            if (abs(p) != 0)  /* Evita operacoes com zero */
                for (l = 1; l <= nv; l++) {
                    if (l != i)
                        Yn[l][j] -= Yn[l][i] * p;
                }
        }
    }
    return 0;

}

int GetNodeNumber(char *, int *, char[MAX_NODES + 1][MAX_NAME + 2]);

void ParametersNewtonRaphsonDiode(double,double,double, double *,double*);

void ParametersNewtonRaphsonEarly(double, double, double,double, double, double, double,
                                  double, double, double *,double*,double*,double*);

int makeNetlistPO(Element[MAX_ELEM], Element[MAX_ELEM],int);

int  addCurrentAboveNode(int *, int *, int ,Element[MAX_ELEM],char[MAX_NODES + 1][MAX_NAME + 2]);

void ListAll(int, int,Element[MAX_ELEM],char[MAX_NODES + 1][MAX_NAME + 2]);

int getLastNFromString(string,int, char*);

void SetMNAPO(int,int,double[MAX_NODES + 1][MAX_NODES + 2], Element[MAX_ELEM],bool,bool,int);

void GetResultsPO(int,Element[MAX_ELEM],double[MAX_NODES+1]);





#endif //UNTITLED_FUNCTIONS_H

void DiodeResistanceandCapacitance(double ,double ,double , double , double , double *, double *);