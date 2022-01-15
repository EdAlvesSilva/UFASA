
/*
Elementos aceitos e linhas do netlist:

Resistor:  R<nome> <no+> <no-> <resistencia>
VCCS:      G<nome> <io+> <io-> <vi+> <vi-> <transcondutancia>
VCVC:      E<nome> <vo+> <vo-> <vi+> <vi-> <ganho de tensao>
CCCS:      F<nome> <io+> <io-> <ii+> <ii-> <ganho de corrente>
CCVS:      H<nome> <vo+> <vo-> <ii+> <ii-> <transresistencia>
Fonte I:   I<nome> <io+> <io-> <corrente>
Fonte V:   V<nome> <vo+> <vo-> <tensao>
Amp. op.:  O<nome> <vo1> <vo2> <vi2> <vi2>

As fontes F e H tem o ramo de entrada em curto
O amplificador operacional ideal tem a saida suspensa
Os nos podem ser nomes
*/

#define versao "1.0j - 21/06/2017"

#include <iostream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <complex>
#include "Functions.h"
#include "psFunctions.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ShowErrorMessage(int error){
    cout <<"Error "<<error<<endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////  MAIN  ////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(void) {

    bool mode,NRmode;
    int error;

    Element netlist[MAX_ELEM]; /* Netlist */
    Element netlistPO[MAX_ELEM];
    Element netlistPS[MAX_ELEM];

    int
            ne, nePO, nePS, /* Elementos */
            nv, nvPO, nvPS, /* Variaveis */
            nn, nnPO, nnPS, /* Nos */
            i, j,
            anal_index,
            indexNR;

    string texto;

    char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
            filename[MAX_LINE + 1],
            type,
            na[MAX_NAME], nb[MAX_NAME], nc[MAX_NAME], nd[MAX_NAME],
            lista[MAX_NODES + 1][MAX_NAME + 2], /*Tem que caber jx antes do nome */
            listaPO[MAX_NODES + 1][MAX_NAME + 2], /*Tem que caber jx antes do nome */
            listaPS[MAX_NODES + 1][MAX_NAME + 2], /*Tem que caber jx antes do nome */
            txt[MAX_LINE + 1],
            anal_type[2],
            *p;


    double
            Yn[MAX_NODES + 1][MAX_NODES + 2],
            YnPO[MAX_NODES + 1][MAX_NODES + 2],
            previousSolution[MAX_NODES+1],
            number_of_frequencies,
            anal_begin,
            anal_end,
            anal_points,
            frequency[MAX_FREQUENCIES],
            value,val_is,val_vt,val_c0,val_c1,module,phase,
            val_alpha,val_alphaR,val_isBe,val_vtBe,val_isBc,val_vtBc,
            val_va,val_c0Be,val_c1Be,val_c0Bc,val_c1Bc;


    ifstream file;
    ofstream savefilePO,savefilePS;
    char * save_filenamePO, * save_filenamePS;

    complex<double> YnPS[MAX_NODES + 1][MAX_NODES + 2];

    number_of_frequencies=10;
    system("cls");
    cout << "Programa de Circuilos Eletricos II\n";
    cout << "Por Carolina Lazzari Bez e Eduardo Alves da Silva\n";
    cout << "Versao" << versao << endl;
    repeat:
    /* Leitura do netlist */
    ne = 0;
    nv = 0;
    strcpy(lista[0], "0");
    cout << "Nome do file com o netlist (ex: mna.net): ";
    cin >> filename;

    cout<<"1\n";
    file.open(filename, ios::in);

    cout<<"2\n";
    try {

        if (not file.is_open()) {
            throw NO_SUCH_FILE_ERROR;
        }
    }
    catch (int error) {
        cout << "Erro " << error << "- Arquivo " << filename << " inexistente.\n";
        goto repeat;

    }






    cout << "Lendo netlist:\n";
    getline(file, texto);
    cout << "Titulo:" << texto << endl;
    for (i=0;i<=number_of_frequencies;i++) {
        frequency[i] = i*10; //inicializa as frequencias a serem analisadas
    }
    for (i=0;i<=MAX_NODES+1;i++){
        for (j=0;j<=MAX_NODES+2;j++){
            Yn[i][j]=0;
            YnPO[i][j]=0;
            //
            // YnPS[i][j]=complex<double> (0,0);
        }
    }

    while (getline(file, texto)) {
        strncpy(txt, texto.c_str(), MAX_LINE);
        ne++; /* Nao usa o netlist[0] */
        if (ne > MAX_ELEM) {
            cout << "O programa so aceita ate" << MAX_ELEM << "elementos\n";
            exit(1);
        }
        txt[0] = char(toupper(txt[0]));
        type = txt[0];
        char name[MAX_NAME];
        char nameLA[MAX_NAME];
        char nameLB[MAX_NAME];
        char qType[MAX_NAME];

        sscanf(txt, "%10s", name);
        netlist[ne].setName(name);
        p = txt + strlen(netlist[ne].getName()); /* Inicio dos parametros */
        /* O que e lido depende do tipo */
        if (type == 'R' ) {
            sscanf(p, "%10s%10s%lg", na, nb, &value);
            netlist[ne].setValue(value);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " << netlist[ne].getValue() << endl;
            netlist[ne].setA(GetNodeNumber(na, &nv, lista));
            netlist[ne].setB(GetNodeNumber(nb, &nv, lista));

        } else if (type == 'I' || type == 'V'){
            sscanf(p, "%10s%10s%lg%lg%lg", na, nb,&module,&phase, &value);
            Source source;
            source.setModule(module);
            source.setPhase(phase);
            source.setValue(value);
            netlist[ne].setA(GetNodeNumber(na, &nv, lista));
            netlist[ne].setB(GetNodeNumber(nb, &nv, lista));
            netlist[ne].setValue(value);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " << source.getModule()<< " "<< source.getPhase()<< " "<< source.getValue() << endl;
            netlist[ne].setSourceElement(&source);
        } else if (type == 'G' || type == 'E' || type == 'F' || type == 'H') {

            sscanf(p, "%10s%10s%10s%10s%lg", na, nb, nc, nd, &value);
            netlist[ne].setA(GetNodeNumber(na,&nv,lista));
            netlist[ne].setB(GetNodeNumber(nb,&nv,lista));
            netlist[ne].setC(GetNodeNumber(nc,&nv,lista));
            netlist[ne].setD(GetNodeNumber(nd,&nv,lista));
            netlist[ne].setValue(value);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " << nc << " " << nd << " "
                 << netlist[ne].getValue() << endl;
        } else if (type == 'O') {
            sscanf(p, "%10s%10s%10s%10s", na, nb, nc, nd);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " << nc << " " << nd << endl;
            netlist[ne].setA(GetNodeNumber(na,&nv,lista));
            netlist[ne].setB(GetNodeNumber(nb,&nv,lista));
            netlist[ne].setC(GetNodeNumber(nc,&nv,lista));
            netlist[ne].setD(GetNodeNumber(nd,&nv,lista));
        } else if (type == 'D') {
            sscanf(p, "%10s%10s%lg%lg%lg%lg", na, nb, &val_is, &val_vt,&val_c0,&val_c1);
            Diode diode;
            netlist[ne].setA(GetNodeNumber(na,&nv,lista));
            netlist[ne].setB(GetNodeNumber(nb,&nv,lista));
            diode.setIs(val_is);
            diode.setVt(val_vt);
            diode.setC0(val_c0);
            diode.setC1(val_c1);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " << diode.getIs() << diode.getVt() << diode.getC0() << diode.getC1() <<" " << nd << endl;
            netlist[ne].setDiodeElement(&diode);
        } else if (type == 'C') {
            sscanf(p, "%10s%10s%lg", na, nb, &value);
            Capacitor capacitor;
            netlist[ne].setA(GetNodeNumber(na,&nv,lista));
            netlist[ne].setB(GetNodeNumber(nb,&nv,lista));
            netlist[ne].setValue(value);
            netlist[ne].setCapacitorElement(&capacitor);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " <<  netlist[ne].getValue() << endl;
        } else if (type == 'L') {
            sscanf(p, "%10s%10s%lg", na, nb, &value);
            Inductor inductor;
            netlist[ne].setA(GetNodeNumber(na, &nv, lista));
            netlist[ne].setB(GetNodeNumber(nb, &nv, lista));
            netlist[ne].setValue(value);
            netlist[ne].setInductorElement(&inductor);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " << netlist[ne].getValue() << endl;
        } else if (type == 'Q'){


            Transistor transistor;
            sscanf(p, "%10s%10s%10s%10s%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg", na, nb, nc,qType,&val_alpha,&val_alphaR,
                   &val_isBe, &val_vtBe,&val_isBc, &val_vtBc,&val_va,&val_c0Be,&val_c1Be,&val_c0Bc,&val_c1Bc);
            if(!val_alpha)
                val_alpha = ALPHA_DEF;
            if(!val_alphaR)
                val_alphaR = ALPHAR_DEF;
            if(!val_isBe)
                val_isBe = IS_DEF;
            if(!val_isBc)
                val_isBc = IS_DEF;
            if(!val_vtBe)
                val_vtBe = VT_DEF;
            if(!val_vtBc)
                val_vtBc = VT_DEF;
            if(!val_va)
                val_va = VA_DEF;
            if(!val_c0Be)
                val_c0Be = CO_DEF;
            if(!val_c0Bc) {
                val_c0Bc = CO_DEF;
            }
            if(!val_c1Be)
                val_c1Be = CI_DEF;
            if(!val_c1Bc)
                val_c1Bc = CI_DEF;
            if (!strcmp(qType,"NPN"))
                transistor.setQType(NPN);
            else if (!strcmp(qType,"PNP"))
                transistor.setQType(PNP);
            else{
                ShowErrorMessage(TRANSISTOR_INVALID_TYPE_ERROR);
                exit(TRANSISTOR_INVALID_TYPE_ERROR);
            }

            netlist[ne].setA(GetNodeNumber(na,&nv,lista));
            netlist[ne].setB(GetNodeNumber(nb,&nv,lista));
            netlist[ne].setC(GetNodeNumber(nc,&nv,lista));
            transistor.setA(GetNodeNumber(na,&nv,lista));
            transistor.setB(GetNodeNumber(nb,&nv,lista));
            transistor.setC(GetNodeNumber(nc,&nv,lista));
            transistor.setIsBe(val_isBe);
            transistor.setVtBe(val_vtBe);
            transistor.setC0Be(val_c0Be);
            transistor.setC1Be(val_c1Be);
            transistor.setAlpha(val_alpha);
            transistor.setAlphaR(val_alphaR);
            transistor.setVa(val_va);
            transistor.setIsBc(val_isBc);
            transistor.setVtBc(val_vtBc);
            transistor.setC0Bc(val_c0Bc);
            transistor.setC1Bc(val_c1Bc);
            cout << netlist[ne].getName() << " " << na << " " << nb << " " << nc << " " << qType << " " <<
                 transistor.getAlpha() << " " << transistor.getAlphaR() << " " <<
                 transistor.getIsBe() << " " << transistor.getVtBe() << " " <<
                 transistor.getIsBc() << " " <<transistor.getVtBc() << " " <<
                 transistor.getVa() << " " <<
                 transistor.getC0Be() << " " <<transistor.getC1Be() <<" " <<
                 transistor.getC0Bc() << " " <<transistor.getC1Bc() <<" " << nd << endl;
            netlist[ne].setTransistorElement(&transistor);

        } else if (type == 'K') {
            Coupling coupling;
            sscanf(p, "%10s%10s%lg", nameLA, nameLB, &value);
            coupling.setNameLA(nameLA);
            coupling.setNameLB(nameLB);
            coupling.setValue(value);
            netlist[ne].setCouplingElement(&coupling);
            cout << netlist[ne].getName() << " " << coupling.getNameLA() << " " << coupling.getNameLB() << " " <<  coupling.getValue() << endl;


        }
        else if (type == '*') { /* Comentario comeca com "*" */
            cout << "Comentario: " << txt;
            ne--;
        }
        else if (type == '.') {
            strcpy(anal_type,strtok(p," "));
            cout << "|" <<  p << endl;
            anal_begin = atof(strtok(NULL," "));
            anal_end = atof(strtok(NULL," "));
            anal_points = atof(strtok(NULL," "));
            if (!(strcmp(anal_type,"LIN"))) {
                anal_index = LINEAR;
                cout << "Tipo de analise : LINEAR" << endl;
            }
            if (!(strcmp(anal_type,"OCT"))) {
                anal_index = LOG_OCT;
                cout << "Tipo de analise : OITAVA" << endl;
            }
            if (!(strcmp(anal_type,"DEC"))) {
                anal_index = LOG_DEC;
                cout << "Tipo de analise : DECADA" << endl;
            }
            ne--;
            cout << "Inicio : " <<anal_begin << " Final : " << anal_end << " PPD : " << anal_points << endl;
        }

        else {
            cout << "Element desconhecido: " << txt << endl;
            std::cin.get();
            exit(1);
        }
    }
    file.close();
    /* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */

    ////////////////////// Modo Ponto de OperaÃ§Ã£o ////////////////////////////////////


    nvPO=nv;
    nePO = makeNetlistPO(netlist,netlistPO,ne);
    for (i = 0; i <= nv; i++) {
        strcpy(listaPO[i], lista[i]);
    }
    error =  addCurrentAboveNode(&nnPO,&nvPO, nePO,netlistPO,listaPO);
    if (error){
        ShowErrorMessage(error);
        exit(error);
    }

    //std::cin.get();
    ListAll(nvPO,nePO,netlistPO,listaPO);
    std::cin.get();
    mode = PO_MODE; /* Comeca pela analise de ponto de operaÃ§Ã£o */
    NRmode = NR_OFF;
    indexNR=0;
    /* Monta o sistema nodal modificado */
    cout << "O circuito de Ponto de Operacao tem " << nnPO << " nos, " << nvPO << " variaveis e " << nePO << " elementos\n";
    std::cin.get();
    solvesystem:
    for (i=0;i<=nvPO;i++) {
        if (NRmode == NR_ON) {
            previousSolution[i] = YnPO[i][nvPO + 1];
        }
    }

    SetMNAPO(nvPO,nePO,YnPO,netlistPO,mode,NRmode,indexNR);


    /* Resolve o sistema */
    if (SolveSystem(nvPO,YnPO)) {
        std::cin.get();
        exit;
    }
    if (NRmode == NR_ON) {
        NRmode=NR_OFF;
        for (i=1;i<=nvPO;i++) {
            if (fabs(previousSolution[i]- YnPO[i][nvPO + 1])>TOLG) {
                NRmode = NR_ON;
                //std::cin.get();
                indexNR++;
                goto solvesystem;
            }
        }

        if (NRmode == NR_OFF)
            goto endNR;
    }
    else {
        NRmode = NR_ON;
        //std::cin.get();
        indexNR++;
        goto solvesystem;
    }
    endNR:
#ifdef DEBUG
    /* Opcional: Mostra o sistema resolvido */
    cout << "Sistema resolvido:\n";
    for (i = 1; i <= nvPO; i++) {
        for (j = 1; j <= nvPO + 1; j++)
            if (YnPO[i][j] != 0) //cout << YnPO[i][j];
                cout << setw(10) << setprecision(5) << showpos << YnPO[i][j];
            else cout << setw(10) << " ... ";
        cout << "\n";
    }
    cin.get();
#endif
    save_filenamePO = strtok(filename, ".");
    strcat(save_filenamePO,"PO.tab");
    savefilePO.open(save_filenamePO, ios::out);
    try {
        if (not savefilePO.is_open()) {
            throw NO_SUCH_FILE_ERROR;
        }
    }
    catch (int error) {
        cout << "Erro " << error << "- Arquivo " << save_filenamePO << " inexistente.\n";

    }
    /* Mostra solucao PO */
    cout << "Solucao Ponto de Operacao:\n";
    strcpy(txt, "Tensao");
    for (i = 1; i <= nvPO; i++) {
        if (i == nnPO + 1) {
            strcpy(txt, "Corrente");
        }
        cout << txt << " " << listaPO[i] << ": " << YnPO[i][nvPO + 1] << endl;
        savefilePO << listaPO[i] << " " << YnPO[i][nvPO + 1] << endl;
    }
    GetResultsPO(ne,netlist,previousSolution);

    savefilePO.close();

    /////////////////////////////// FIM MODO PONTO DE OPERACAO ////////////////////////////////

    nvPS = nv;
    for (i = 0; i <= nv; i++) {
        strcpy(listaPS[i], lista[i]);
    }
    nePS = makeNetlistPS(netlist,netlistPS,ne,previousSolution);
    error =  addCurrentAboveNode(&nnPS,&nvPS, nePS,netlistPS,listaPS);
    if (error){
        ShowErrorMessage(error);
        exit(error);
    }

    ListAll(nvPS,nePS,netlistPS,listaPS);
    std::cin.get();

    cout << "O circuito de Pequenos Sinais tem " << nnPS << " nos, " << nvPS << " variaveis e " << nePS << " elementos\n";

    for (i=1;i<=nePS;i++) {
        cout << netlistPS[i].getName() << "///" << netlistPS[i].getValue() << endl;
    }
    //////// Calculo das frequencias ////////
    makeFrequenciesList(anal_begin,anal_end,anal_points,frequency,&number_of_frequencies,anal_index);

    int l;
    for (i=0;i<=number_of_frequencies;i++) {
        cout << frequency[i] << " , " << i << endl;
    }
    cin.get();


    cout << filename << endl << endl << endl << endl <<endl;

    save_filenamePS = strtok(filename, ".");
    strcat(save_filenamePS,"PS.tab");
    savefilePS.open(save_filenamePS, ios::out);
    try {
        if (not savefilePS.is_open()) {
            throw NO_SUCH_FILE_ERROR;
        }
    }
    catch (int error) {
        cout << "Erro " << error << "- Arquivo " << save_filenamePS << " inexistente.\n";

    }

    savefilePS << "f\t";
    for (i=1;i<=nvPS;i++) {
        savefilePS << listaPS[i] << "m\t" << listaPS[i] << "f\t";
    }
    savefilePS << endl;
    i = 0;
    for (i=0;i<=number_of_frequencies;i++) {
        cout << "Frequencia analisada :" << frequency[i] << endl;
        SetMNAPS(nvPS,nePS,YnPS,netlistPS,frequency[i]);
        error = SolveSystem(nvPS,YnPS);
        if (error){
            ShowErrorMessage(error);
            exit(error);
        }
        int k;
        savefilePS << frequency[i] << "\t";
        for (k=1;k<=nvPS;k++) {
            for (l=1;l<=nvPS+1;l++) {
                cout <<"[" << scientific <<setprecision(3) << YnPS[k][l].real() << "+ j" << YnPS[k][l].imag() << "]";
            }
            cout << endl;
        }
        /* Mostra solucao PS */
        cout << "Solucao PS:\n";
        strcpy(txt, "Tensao");

        for (k = 1; k <= nvPS; k++) {
            if (k >= nnPS ) {
                strcpy(txt, "Corrente");
            }
            cout << txt << " " << listaPS[k] << ": " << YnPS[k][nvPS+1].real() << "+ j" << YnPS[k][nvPS+1].imag() << endl;
            savefilePS << setprecision(3) << abs(YnPS[k][nvPS+1]) << "\t" << (arg(YnPS[k][nvPS+1])*180)/PI << "\t";
        }
        savefilePS << endl;
    }

    savefilePS.close();

    cout<< "Aperte qualquer tecla para analisar outro circuito. Tecla q ou ESC para sair.\n";
    int c;
    c = cin.get();
    if (not(c=='q'||c==27))
        goto repeat;
    return 0;
}

