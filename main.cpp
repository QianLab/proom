//
//  main.cpp
//  proom
//
//  Pessimistic ROOM (Relaxed LP version for ROOM)
//  Created by Meltem Apaydin on 1/11/16.
//  Copyright (c) 2016 Meltem Apaydin. All rights reserved.
//

#include <stdio.h>
//#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ilcplex/ilocplex.h>

using namespace std;

struct subject { // For Stoichiometric matrix
    vector<string> metabolite;
    vector<string> reaction;
    vector<double> stoic;
    vector<int> mapped_number_index;
    vector<int> mapped_number_index_cmp;
    
};

struct parameter { // Data Parameters
    vector<int> metab_name_index;
    vector<string> metabolite_names;
    vector<string> reaction_names;
    vector<int> reaction_types;
    vector<double> antcore_max;
    vector<double> antcore_min;
    vector<double> wild_type;
    double M = 1000;
    int K = 3; // Allaowable Knockout Number, will be defined by user
    double minbiomass = 5; // If minbiomass change, UPDATE MIN, MAX AND w(j) VALUES FOR REAC FLUX. !!!!!
    double ebs = 0.1; // Defined by user
    double glc_uptake = 100 ;// Glucose uptake rate
};

ILOSTLBEGIN //ILOG standard template library

int main(int argc, const char * argv[]) {
    
    ifstream inFile;
    parameter Data;
    subject S_ij;
    
    S_ij.mapped_number_index_cmp.resize(400);
    
    /////////////////////////////// SET OF METABOLITES, NAMES /////////////////////////////////////
    inFile.open("/Users/meltemapaydin/Desktop/cplex00/cplex00/AntCore_cmp.txt"); // open file
    if (inFile.good()) {

        string line;

        while (!inFile.eof()) {
            getline(inFile,line);
            Data.metabolite_names.push_back(line);
        }
    }
    inFile.close();
    
    ///////////////////////////////// SET OF REACTIONS, NAMES //////////////////////////////////////
    inFile.open("/Users/meltemapaydin/Desktop/cplex00/cplex00/AntCore_rxnnames.txt");
    if (inFile.good()) {
        
        string line;
        while (!inFile.eof()) {
            getline(inFile,line);
            Data.reaction_names.push_back(line);
        }
//        vector<string>::iterator it;
//        for (it = rxn_namesVector.begin(); it < rxn_namesVector.end(); it++) {
//            cout << *it << endl;
//        }
    }
    inFile.close();

    ////////////////////////////////// REACTION TYPES (integer vals) ////////////////////////////////
    inFile.open("/Users/meltemapaydin/Desktop/cplex00/cplex00/AntCore_rxntype.txt");
    if (inFile.good()) {
        
        string line;
        while (!inFile.eof()) {
            getline(inFile,line);
            int reac_type;
            string reac_name;
            istringstream iss(line);
            iss >> reac_name >> reac_type;
            Data.reaction_types.push_back(reac_type);
        }
        //for (const auto& p : rxn_typeVector)
        //{
        //    cout << p.first << "  -----> " << p.second << std::endl;
        //    // or std::cout << std::get<0>(p) << ", " << std::get<1>(p) << std::endl;
        //}
    }
    inFile.close();


    //////////////////////////////////////// S_ij ////////////////////////////////////////////
    inFile.open("/Users/meltemapaydin/Desktop/cplex00/cplex00/AntCore_sij.txt"); // open file
    if (inFile.good()) {
        
        string line;
        while (!inFile.eof()) {
            getline(inFile,line);
            string reac_name;
            double val;
            size_t pos0 = line.find(".");
            string str0 = line.substr(0,pos0);
            string str1 = line.substr(pos0+1);
            istringstream iss(str1);
            
            iss >> reac_name >> val ;
            
            S_ij.metabolite.push_back(str0);
            S_ij.reaction.push_back(reac_name);
            S_ij.stoic.push_back(val);
            
        }
        S_ij.metabolite.erase (S_ij.metabolite.begin());
        S_ij.metabolite.erase(S_ij.metabolite.end()-1);
        S_ij.metabolite.pop_back();
        
        S_ij.reaction.erase (S_ij.reaction.begin());
        S_ij.reaction.erase(S_ij.reaction.end()-1);
        S_ij.reaction.pop_back();
        
        S_ij.stoic.erase (S_ij.stoic.begin());
        S_ij.stoic.erase(S_ij.stoic.end()-1);
        S_ij.stoic.pop_back();
        
        int sayac = 0;
        S_ij.mapped_number_index.push_back(sayac);
        
        for (vector <string>::iterator ff = S_ij.reaction.begin() ; ff < S_ij.reaction.end()-1 ; ff++){
            
            if (*ff == *(ff+1)) {
                S_ij.mapped_number_index.push_back(sayac);
            }
            else if (*ff != *(ff+1)){
                sayac++;
                S_ij.mapped_number_index.push_back(sayac);
            }
        }
        //        for(int i=0; i<S_ij.metabolite.size(); i++)
        //        {
        //          cout << S_ij.metabolite[i] << endl;
        //        }
    }
    inFile.close();
    
    ////////////////////////////////////// AntCoreMAx /////////////////////////////////////////
    inFile.open("/Users/meltemapaydin/Desktop/cplex00/cplex00/AntCoreMax.txt");
    if (inFile.good()) {
        
        
        string line;
        while (!inFile.eof()) {
            getline(inFile,line);
            istringstream iss(line);
            string p1;
            double p2;
            iss >> p1 >> p2;
            Data.antcore_max.push_back(p2);
        }
    }
    inFile.close();
    
    ////////////////////////////////////////// AntCoreMin //////////////////////////////////////
    inFile.open("/Users/meltemapaydin/Desktop/cplex00/cplex00/AntCoreMin.txt");
    if (inFile.good()) {
     
        string line;
        while (!inFile.eof()) {
            getline(inFile,line);
            istringstream iss(line);
            string p1;
            double p2;
            iss >> p1 >> p2;
            Data.antcore_min.push_back(p2);
        }
    }
    inFile.close();
    
    //////////////////////////////////// Wild Type //////////////////////////////////////////////
    inFile.open("/Users/meltemapaydin/Desktop/cplex00/cplex00/w_j.txt");
    if (inFile.good()) {
 
        string line;
        while (!inFile.eof()) {
            getline(inFile,line);
            istringstream iss(line);
            string p1;
            double p2;
            iss >> p1 >> p2;
            Data.wild_type.push_back(p2);
        }
    }
    inFile.close();
    
    Data.reaction_names.erase (Data.reaction_names.begin());
    Data.reaction_names.erase(Data.reaction_names.end()-1);
    Data.reaction_names.pop_back();
    
    Data.reaction_types.erase (Data.reaction_types.begin());
    Data.reaction_types.erase(Data.reaction_types.end()-1);
    Data.reaction_types.pop_back();
    
    Data.antcore_max.erase (Data.antcore_max.begin());
    Data.antcore_max.erase(Data.antcore_max.end()-1);
    Data.antcore_max.pop_back();
    
    Data.antcore_min.erase (Data.antcore_min.begin());
    Data.antcore_min.erase(Data.antcore_min.end()-1);
    Data.antcore_min.pop_back();
    
    Data.wild_type.erase (Data.wild_type.begin());
    Data.wild_type.erase(Data.wild_type.end()-1);
    Data.wild_type.pop_back();
    
    Data.metabolite_names.erase (Data.metabolite_names.begin());
    Data.metabolite_names.erase(Data.metabolite_names.end()-1);
    Data.metabolite_names.pop_back();
    
    
    //////////////// CPLEX PART //////////////////////
    IloEnv env;
    try {

    IloModel model(env);
    
    IloInt nReac = Data.reaction_names.size();
    IloInt nMetab = Data.metabolite_names.size();
    
    for (int sy = 0 ; sy<nMetab ; sy++){
        Data.metab_name_index.push_back(sy);
    }
    
    for (vector<int>::iterator it = Data.metab_name_index.begin(); it < Data.metab_name_index.end(); it++){
        string temp = Data.metabolite_names[it-Data.metab_name_index.begin()];
        for (vector<string>::iterator it1 = S_ij.metabolite.begin(); it1<S_ij.metabolite.end();it1++ ){
            size_t foundd = temp.find(*it1,0);
            if (foundd != string::npos){
                S_ij.mapped_number_index_cmp[it1-S_ij.metabolite.begin()]= *it;
            }
        }
    }

    
//    cout << "number of metabolites: " << nMetab << endl;
//    cout << "Give knockout number: " ;
//    cin >> Data.K;
//    
//    cout << "Give Epsilon value between 0-1: " ;
//    cin >> Data.ebs;
    
    // CONTINUOUS DEC. VARIABLES
    
    char VarName[24];
    
    IloNumVarArray flux_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName, "v(%ld)", j);
        //flux_values.add(IloNumVar(env, Data.antcore_min[j], Data.antcore_max[j] , ILOFLOAT , VarName));
        flux_values.add(IloNumVar(env , VarName));
    }
    
    char VarNameup[24];
    
    IloNumVarArray upflux_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNameup, "v_up(%ld)", j);
        upflux_values.add(IloNumVar(env , VarNameup));
    }
    
    IloNumVar uglc(env, "uglc");
    
    char VarNamecmp[24];
    
    IloNumVarArray dual_metabolite(env); // u(i): dual variables for stoich. constraint
    for (IloInt i = 0; i < nMetab; i++){
        sprintf( VarNamecmp, "u(%ld)", i);
        dual_metabolite.add(IloNumVar(env , VarNamecmp));
    }
    
    IloNumVar ubiom(env, 0.0, IloInfinity, "ubiom"); // nonneg. dual var for biomass constraint
    
    char VarNamemin[24];
    
    IloNumVarArray umin_values(env); // nonneg. dual variables for lower bound constraints
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNamemin, "umin(%ld)", j);
        umin_values.add(IloNumVar(env , 0.0, IloInfinity, VarNamemin));
    }
    
    char VarNamemax[24];
    
    IloNumVarArray umax_values(env); // nonneg. dual variables for upper bound constraints
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNamemax, "umax(%ld)", j);
        umax_values.add(IloNumVar(env , 0.0, IloInfinity, VarNamemax));
    }
    ////
    char VarNamemin2[24];
    
    IloNumVarArray umin2_values(env); // nonneg. dual variables for ROOM ineq. constraints: umin2(j)
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNamemin2, "umin2(%ld)", j);
        umin2_values.add(IloNumVar(env , 0.0, IloInfinity, VarNamemin2));
    }
    
    char VarNamemax2[24];
    
    IloNumVarArray umax2_values(env); // nonneg. dual variables for ROOM ineq. constraints: umax2(j)
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNamemax2, "umax2(%ld)", j);
        umax2_values.add(IloNumVar(env , 0.0, IloInfinity, VarNamemax2));
    }
    //
    char VarNamea[24];
    
    IloNumVarArray a_values(env); // nonneg. dual variables for 0<y(j): a(j)
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNamea, "a(%ld)", j);
        a_values.add(IloNumVar(env , 0.0, IloInfinity, VarNamea));
    }
    
    char VarNameb[24];
    
    IloNumVarArray b_values(env); // nonneg. dual variables for y(j)<1: b(j)
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNameb, "b(%ld)", j);
        b_values.add(IloNumVar(env , 0.0, IloInfinity, VarNameb));
    }
    
    IloNumVar c(env, 0.0, IloInfinity, "c"); // nonneg. dual var for relaxation constr : c
    
    char VarNamey[24];
    
    IloNumVarArray y_values(env); // y(j) between 0 and 1
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNamey, "y(%ld)", j);
        y_values.add(IloNumVar(env , 0.0, 1.0, VarNamey));
    }
    
    char VarNamey_up[24];
    
    IloNumVarArray y_up_values(env); // y_up(j) between 0 and 1
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarNamey_up, "y_up(%ld)", j);
        y_up_values.add(IloNumVar(env , 0.0, 1.0, VarNamey_up));
    }
    
    //////// BINARY DEC VARIABLES /////
    
    char VarName_knockout[24];
    
    IloBoolVarArray z_values(env); // z(j): knockout binary  dec. variable
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName_knockout, "z(%ld)", j);
        z_values.add(IloBoolVar(env , VarName_knockout));
    }
    
    IloBoolVar bin1(env, "bin1");
    IloBoolVar bin7(env, "bin7");
    
    char VarName_bin2[24];
    IloBoolVarArray bin2_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName_bin2, "bin2(%ld)", j);
        bin2_values.add(IloBoolVar(env , VarName_bin2));
    }
    
    char VarName_bin3[24];
    IloBoolVarArray bin3_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName_bin3, "bin3(%ld)", j);
        bin3_values.add(IloBoolVar(env , VarName_bin3));
    }
    
    char VarName_bin4[24];
    IloBoolVarArray bin4_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName_bin4, "bin4(%ld)", j);
        bin4_values.add(IloBoolVar(env , VarName_bin4));
    }
    
    char VarName_bin5[24];
    IloBoolVarArray bin5_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName_bin5, "bin5(%ld)", j);
        bin5_values.add(IloBoolVar(env , VarName_bin5));
    }
    
    char VarName_bin6[24];
    IloBoolVarArray bin6_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName_bin6, "bin6(%ld)", j);
        bin6_values.add(IloBoolVar(env , VarName_bin6));
    }
    
    char VarName_bin8[24];
    IloBoolVarArray bin8_values(env);
    for (IloInt j = 0; j < nReac; j++){
        sprintf( VarName_bin8, "bin8(%ld)", j);
        bin8_values.add(IloBoolVar(env , VarName_bin8));
    }
    
    //////// OBJECTIVE FUNCTION ////////
    
    // build objective function expression
    IloExpr exprObj(env);
    exprObj = flux_values[110];
    // add obj function to model
    model.add(IloMaximize(env, exprObj));
    exprObj.end();
    
    ///////////////////////////// EXPRESSIONS FOR CONSTRAINTS ///////////////////////////
    
    ////////// Updating the LB and UB according to reaction types //////////
    vector<int> UB(nReac);
    vector<int> LB(nReac);
    
    for (vector<int>::iterator i = Data.reaction_types.begin() ; i < Data.reaction_types.end() ; i++ ){
        
        if (*i == 0){
            LB[i-Data.reaction_types.begin()] = 0;
            UB[i-Data.reaction_types.begin()] = 1000;
        }
        else if (*i == 1){
            LB[i-Data.reaction_types.begin()] = -1000;
            UB[i-Data.reaction_types.begin()] = 1000;
        }
        else if (*i == 2){
            LB[i-Data.reaction_types.begin()] = 0;
            UB[i-Data.reaction_types.begin()] = 0;
        }
        else if (*i == 4){
            LB[i-Data.reaction_types.begin()] = 0;
            UB[i-Data.reaction_types.begin()] = 1000;
        }
        else if (*i == 3){
            LB[i-Data.reaction_types.begin()] = -1000;
            UB[i-Data.reaction_types.begin()] = 1000;
        }
        else cout << "Error in rxn types for reac " << i-Data.reaction_types.begin() << endl;
        
        //cout << "UB " << i-Data.reaction_types.begin() << " for reac_type " << *i << ": " << UB[i-Data.reaction_types.begin()] << endl;
        //cout << "LB " << i-Data.reaction_types.begin() << " for reac_type " << *i << ": " << LB[i-Data.reaction_types.begin()] << endl;
    }
    string temp;
    int temp_pos;
    temp = "'EX_gluc'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    LB[temp_pos] = -100;
    temp = "'EX_o2'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    LB[temp_pos] = -100;
    temp = "'EX_so4'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    LB[temp_pos] = -100;
    temp = "'EX_nh3'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    LB[temp_pos] = -100;
    temp = "'EX_cit'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    LB[temp_pos] = -100;
    temp = "'EX_glyc'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    LB[temp_pos] = -100;
    
    //cout << "size of UB vector: " << UB.size() << "  size of LB:  " << LB.size() << endl;

    /// set LB and UB for flux values and upflux values //////
    for (int i=0 ; i<UB.size() ; i++){
        flux_values[i].setBounds(LB[i], UB[i]);
        upflux_values[i].setBounds(LB[i], UB[i]);
    }
    
    vector<int> blocked(nReac);
    
    for (int i=0; i<LB.size(); i++){
        if (UB[i]==0 && LB[i]==0) {
            blocked[i]=0;
            z_values[i].setBounds(1, 1);
            cout << "i_0: " << i << endl;
            for (vector<int>::iterator it=S_ij.mapped_number_index.begin() ; it<S_ij.mapped_number_index.end() ; it++){
                if (*it==i){
                    //cout << "Position: " << it-S_ij.mapped_number_index.begin() << endl;
                    S_ij.stoic[it-S_ij.mapped_number_index.begin()] = 0;
                    //cout << "stoic val: " << S_ij.stoic[it-S_ij.mapped_number_index.begin()] << endl;
                }
            }
        }
        
        if (UB[i]>0 && LB[i]>0){
            cout << "i_gt: " << i << endl;
            z_values[i].setBounds(1, 1);
        }
        if (UB[i]<0 && LB[i]<0){
            cout << "i_lt: " << i << endl;
            z_values[i].setBounds(1, 1);
        }
    }
    // PRINT OUT PARAMETER VALUES
    //#if 0
    for(vector<double>::iterator i=S_ij.stoic.begin(); i<S_ij.stoic.end(); i++)
    {
        cout << *i << endl;
    }
    //#endif
    
    //////// Outer Problem Constraints ////////////
    
    IloExpr allowable_knockout(env); // Allowable knockout constraint sum( j , (1-z(j)) ) <= K
    for (IloInt j = 0 ; j < nReac ; j++ ){
        allowable_knockout += ( 1 - z_values[j] );
    }
    model.add(allowable_knockout <= Data.K);
    allowable_knockout.end();
    
    ///////////
    
    vector<string>::iterator it;
    string::size_type s;
    
    char stoic_cnstrt_name[24];
    IloRangeArray stoi_cst(env,nMetab);
    int count = 0;
    
    for (it = Data.metabolite_names.begin() ; it < Data.metabolite_names.end(); it++){
        
        IloExpr expr_stoic(env);
        string temp = *it;
        int notfind = 0;
        for (vector<string>::iterator it1=S_ij.metabolite.begin(); it1< S_ij.metabolite.end(); it1++){
            
            s = temp.find(*it1,0);
            
            if( s != string::npos ){
                notfind++;
                //auto pos = it1 - S_ij.metabolite.begin();
                //cout << "Pos: " << pos <<"  SONUC: "<< S_ij.metabolite[it1 - S_ij.metabolite.begin()] << " stoic value: " << S_ij.stoic[it1 - S_ij.metabolite.begin()] << "  corres v: " <<S_ij.mapped_number_index[it1 - S_ij.metabolite.begin()] <<endl ;
                expr_stoic += S_ij.stoic[it1 - S_ij.metabolite.begin()] * upflux_values[(IloInt)S_ij.mapped_number_index[it1 - S_ij.metabolite.begin()]];
            }
        }
        sprintf(stoic_cnstrt_name, "stoi_cst(%d)",count+1);
        
        if (notfind != 0){
            count++;
            stoi_cst.add(IloRange(env, 0 , expr_stoic , 0 , stoic_cnstrt_name));
        }
        expr_stoic.end();
    }
    model.add(stoi_cst);
    stoi_cst.end();
    
    //////
    temp = "'EX_gluc'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    IloExpr glucose(env); // fixed glucose rate to 100
    glucose  = upflux_values[temp_pos] ;
    model.add(glucose == -100);
    glucose.end();
    
    ///////
    temp = "'75'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    IloExpr biomass(env); // fixed glucose rate to 100
    biomass  = upflux_values[temp_pos] ;
    model.add(biomass >= Data.minbiomass);
    biomass.end();
    
    ////// LB ve UB leri AYRI EKLEMEYI unutma, ayrica reac_type 2 ise z(j)=1 i de yap
    IloRangeArray upbounds_max(env,nReac);
    char temp_name[24];
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = upflux_values[i] - Data.antcore_max[i] * z_values[i] ;
        sprintf(temp_name,"upbounds_max(%d)",int(i+1));
        upbounds_max.add(IloRange(env,-IloInfinity,xpr,0,temp_name));
        xpr.end();
    }
    model.add(upbounds_max);
    upbounds_max.end();
    
    //////
    IloRangeArray upbounds_min(env,nReac);
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = - upflux_values[i] + Data.antcore_min[i] * z_values[i] ;
        sprintf(temp_name,"upbounds_min(%d)",int(i+1));
        upbounds_min.add(IloRange(env,-IloInfinity,xpr,0,temp_name));
        xpr.end();
    }
    model.add(upbounds_min);
    upbounds_min.end();
    
    ////// ROOM ineq constraiint (distance constraint replicated outer level)
    IloRangeArray roomdist_max(env,nReac);
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = upflux_values[i] - y_up_values[i] * ( Data.antcore_max[i] - Data.wild_type [i] ) ;
        sprintf(temp_name,"dist_max(%d)",int(i+1));
        roomdist_max.add(IloRange(env,-IloInfinity,xpr,Data.wild_type[i],temp_name));
        xpr.end();
    }
    model.add(roomdist_max);
    roomdist_max.end();
    
    //////////
    IloRangeArray roomdist_min(env,nReac);
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = upflux_values[i] - y_up_values[i] * ( Data.antcore_min[i] - Data.wild_type [i] ) ;
        sprintf(temp_name,"dist_min(%d)",int(i+1));
        roomdist_min.add(IloRange(env,Data.wild_type[i],xpr,IloInfinity,temp_name));
        xpr.end();
    }
    model.add(roomdist_min);
    roomdist_min.end();
    
    //////////////////// Constraints for Primal Inner Problem ///////////////////////////
    
    vector<string>::iterator itt;
    string::size_type ss;
    
    //char stoic_cnstrt_name[24];
    IloRangeArray stoi_cst_primal(env,nMetab);
    count = 0;
    
    for (itt = Data.metabolite_names.begin() ; itt < Data.metabolite_names.end(); itt++){
        
        IloExpr expr_stoic(env);
        string temp = *itt;
        int notfind = 0;
        for (vector<string>::iterator it1=S_ij.metabolite.begin(); it1< S_ij.metabolite.end(); it1++){
            
            ss = temp.find(*it1,0);
            
            if( ss != string::npos ){
                notfind++;
                //auto pos = it1 - S_ij.metabolite.begin();
                //cout << "Pos: " << pos <<"  SONUC: "<< S_ij.metabolite[it1 - S_ij.metabolite.begin()] << " stoic value: " << S_ij.stoic[it1 - S_ij.metabolite.begin()] << "  corres v: " <<S_ij.mapped_number_index[it1 - S_ij.metabolite.begin()] <<endl ;
                expr_stoic += S_ij.stoic[it1 - S_ij.metabolite.begin()] * flux_values[(IloInt)S_ij.mapped_number_index[it1 - S_ij.metabolite.begin()]];
            }
        }
        sprintf(stoic_cnstrt_name, "stoi_cst_inner(%d)",count+1);
        
        if (notfind != 0){
            count++;
            stoi_cst_primal.add(IloRange(env, 0 , expr_stoic , 0 , stoic_cnstrt_name));
        }
        expr_stoic.end();
    }
    model.add(stoi_cst_primal);
    stoi_cst_primal.end();

    /////////

    temp = "'EX_gluc'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    IloExpr glucose_inner(env); // fixed glucose rate to 100
    glucose_inner  = flux_values[temp_pos] ;
    model.add(glucose_inner == -100);
    glucose_inner.end();
    
    ////
    temp = "'75'";
    temp_pos =  S_ij.mapped_number_index[find(S_ij.reaction.begin(),S_ij.reaction.end(),temp)-S_ij.reaction.begin()];
    IloExpr biomass_inner(env); // fixed glucose rate to 100
    biomass_inner  = flux_values[temp_pos] ;
    model.add(biomass_inner >= Data.minbiomass);
    biomass_inner.end();
    
    /////
    IloRangeArray bounds_max(env,nReac);
    char max_name[24];
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = flux_values[i] - Data.antcore_max[i] * z_values[i] ;
        sprintf(max_name,"bounds_max(%d)",int(i+1));
        bounds_max.add(IloRange(env,-IloInfinity,xpr,0,max_name));
        xpr.end();
    }
    model.add(bounds_max);
    bounds_max.end();
    
    /////
    IloRangeArray bounds_min(env,nReac);
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = - flux_values[i] + Data.antcore_min[i] * z_values[i] ;
        sprintf(temp_name,"bounds_min(%d)",int(i+1));
        bounds_min.add(IloRange(env,-IloInfinity,xpr,0,temp_name));
        xpr.end();
    }
    model.add(bounds_min);
    bounds_min.end();
    
    //////
    IloRangeArray inner_roomdist_max(env,nReac);
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = flux_values[i] - y_values[i] * ( Data.antcore_max[i] - Data.wild_type [i] ) ;
        sprintf(temp_name,"dist_max_inner(%d)",int(i+1));
        inner_roomdist_max.add(IloRange(env,-IloInfinity,xpr,Data.wild_type[i],temp_name));
        xpr.end();
    }
    model.add(inner_roomdist_max);
    inner_roomdist_max.end();
    
    /////////
    
    IloRangeArray inner_roomdist_min(env,nReac);
    for (IloInt i = 0; i < Data.reaction_names.size(); i++ ){
        IloExpr xpr(env);
        xpr = flux_values[i] - y_values[i] * ( Data.antcore_min[i] - Data.wild_type [i] ) ;
        sprintf(temp_name,"dist_min_inner(%d)",int(i+1));
        inner_roomdist_min.add(IloRange(env,Data.wild_type[i],xpr,IloInfinity,temp_name));
        xpr.end();
    }
    model.add(inner_roomdist_min);
    inner_roomdist_min.end();
    
    ///////
    
    IloExpr relaxation(env); // relaxation constraint
    IloExpr relaxation1(env);
    
    for (IloInt j = 0 ; j < nReac ; j++ ){
        relaxation += y_values[j];
        relaxation1 += y_up_values[j];
    }
    relaxation = relaxation - (1+Data.ebs)*relaxation1;
    model.add(relaxation <= 0);
    relaxation1.end();
    relaxation.end();
    
    ///////////////////////////////// KKT CONDITIONS OF THE INNER PROBLEM //////////////////////////////////
    IloRangeArray kkt1(env,nReac);
    IloRangeArray kkt2(env,nReac);
    IloRangeArray kkt3(env,nReac);
    IloRangeArray kkt4(env,nReac);
    IloRangeArray kkt5(env,nReac);
    IloRangeArray kkt6(env,nReac);
    IloRangeArray kkt9(env,nReac);
    IloRangeArray kkt10(env,nReac);
    IloRangeArray kkt11(env,nReac);
    IloRangeArray kkt12(env,nReac);
    IloRangeArray kkt13(env,nReac);
    IloRangeArray kkt14(env,nReac);
    IloRangeArray kkt15(env,nReac);
    IloRangeArray kkt16(env,nReac);
    
    char namme [24];
    
    for (IloInt j = 0 ; j < Data.reaction_names.size() ; j++){
        IloExpr k1(env);
        IloExpr k2(env);
        IloExpr k3(env);
        IloExpr k4(env);
        IloExpr k5(env);
        IloExpr k6(env);
        IloExpr k9(env);
        IloExpr k10(env);
        IloExpr k11(env);
        IloExpr k12(env);
        IloExpr k13(env);
        IloExpr k14(env);
        IloExpr k15(env);
        IloExpr k16(env);
        
        if (j==98){ // v('75')
            temp = Data.reaction_names[98];
            IloExpr subexpr0(env);
            for (vector<string>::iterator it=S_ij.reaction.begin(); it< S_ij.reaction.end(); it++){
                ss = temp.find(*it,0);
                if( ss != string::npos ){
                    int s1 = S_ij.mapped_number_index_cmp[it-S_ij.reaction.begin()];
                    subexpr0 += S_ij.stoic[it-S_ij.reaction.begin()] * dual_metabolite[s1];
                }
            }
            k2 = -umax2_values[j]*(Data.antcore_max[j]-Data.wild_type[j]) + umin2_values[j]*(Data.antcore_min[j]-Data.wild_type[j]) - a_values[j] + b_values[j] + c ;
            k1 = -ubiom - umin_values[j] + umax_values[j] + umax2_values[j] -umin2_values[j]+subexpr0;
        }
        
        else if (j==103){ // v('EX_gluc')
            temp = Data.reaction_names[103];
            IloExpr subexpr1(env);
            for (vector<string>::iterator it=S_ij.reaction.begin(); it< S_ij.reaction.end(); it++){
                ss = temp.find(*it,0);
                if( ss != string::npos ){
                    int s1 = S_ij.mapped_number_index_cmp[it-S_ij.reaction.begin()];
                    subexpr1 += S_ij.stoic[it-S_ij.reaction.begin()] * dual_metabolite[s1];
                }
            }
            k2 = -umax2_values[j]*(Data.antcore_max[j]-Data.wild_type[j]) + umin2_values[j]*(Data.antcore_min[j]-Data.wild_type[j]) - a_values[j] + b_values[j] + c ;
            k1 = uglc - umin_values[j] + umax_values[j] + umax2_values[j] -umin2_values[j]+subexpr1;
        }
        
        else if (j==110){ // v('EX_suc')
            temp = Data.reaction_names[110];
            IloExpr subexpr2(env);
            for (vector<string>::iterator it=S_ij.reaction.begin(); it< S_ij.reaction.end(); it++){
                ss = temp.find(*it,0);
                if( ss != string::npos ){
                    int s1 = S_ij.mapped_number_index_cmp[it-S_ij.reaction.begin()];
                    subexpr2 += S_ij.stoic[it-S_ij.reaction.begin()] * dual_metabolite[s1];
                }
            }
            k2 = -umax2_values[j]*(Data.antcore_max[j]-Data.wild_type[j]) + umin2_values[j]*(Data.antcore_min[j]-Data.wild_type[j]) - a_values[j] + b_values[j] + c ;
            k1 = 1 - umin_values[j] + umax_values[j] + umax2_values[j] -umin2_values[j]+subexpr2;
        }
        else {
            temp = Data.reaction_names[j];
            IloExpr subexpr(env);
            for (vector<string>::iterator it=S_ij.reaction.begin(); it< S_ij.reaction.end(); it++){
                ss = temp.find(*it,0);
                if( ss != string::npos ){
                    IloInt s1 = S_ij.mapped_number_index_cmp[it-S_ij.reaction.begin()];
                    subexpr += S_ij.stoic[it-S_ij.reaction.begin()] * dual_metabolite[s1];
                }
                //else subexpr = 0;
            }
            k2 = -umax2_values[j]*(Data.antcore_max[j]-Data.wild_type[j]) + umin2_values[j]*(Data.antcore_min[j]-Data.wild_type[j]) - a_values[j] + b_values[j] + c ;
            k1 = - umin_values[j] + umax_values[j] + umax2_values[j] -umin2_values[j] + subexpr;
            
        }
        sprintf(namme,"kkt2(%d)",int(j));
        kkt2.add(IloRange(env,0,k2,0,namme));
        k2.end();
        sprintf(namme,"kkt1(%d)",int(j));
        kkt1.add(IloRange(env,0,k1,0,namme));
        k1.end();
        
        /// adding constraints from kkt3 to kkt18
        k3 = umin_values[j] - Data.M*bin2_values[j];
        sprintf(namme,"kkt3(%d)",int(j));
        kkt3.add(IloRange(env,-IloInfinity,k3,0,namme));
        k3.end();
        
        k4 = flux_values[j] - Data.antcore_min[j]*z_values[j] - Data.M*(1-bin2_values[j]);
        sprintf(namme,"kkt4(%d)",int(j));
        kkt4.add(IloRange(env,-IloInfinity,k4,0,namme));
        k4.end();
        
        k5 = umax_values[j] - Data.M*bin3_values[j];
        sprintf(namme,"kkt5(%d)",int(j));
        kkt5.add(IloRange(env,-IloInfinity,k5,0,namme));
        k5.end();
        
        k6 = Data.antcore_max[j]*z_values[j] - flux_values[j]- Data.M*(1-bin3_values[j]);
        sprintf(namme,"kkt6(%d)",int(j));
        kkt6.add(IloRange(env,-IloInfinity,k6,0,namme));
        k6.end();
        
        k9 = umax2_values[j] - Data.M * bin4_values[j];
        sprintf(namme,"kkt9(%d)",int(j));
        kkt9.add(IloRange(env,-IloInfinity,k9,0,namme));
        k9.end();
        
        k10 = Data.wild_type[j] - flux_values[j] + y_values[j]*(Data.antcore_max[j]-Data.wild_type[j]) - Data.M*(1-bin4_values[j]);
        sprintf(namme,"kkt10(%d)",int(j));
        kkt10.add(IloRange(env,-IloInfinity,k10,0,namme));
        k10.end();
        
        k11 = umin2_values[j] - Data.M * bin5_values[j];
        sprintf(namme,"kkt11(%d)",int(j));
        kkt11.add(IloRange(env,-IloInfinity,k11,0,namme));
        k11.end();
        
        k12 = flux_values[j] - y_values[j]*(Data.antcore_min[j]-Data.wild_type[j]) - Data.wild_type[j] - Data.M*(1-bin5_values[j]);
        sprintf(namme,"kkt12(%d)",int(j));
        kkt12.add(IloRange(env,-IloInfinity,k12,0,namme));
        k12.end();
        
        k13 = a_values[j] - Data.M * bin8_values[j];
        sprintf(namme,"kkt13(%d)",int(j));
        kkt13.add(IloRange(env,-IloInfinity,k13,0,namme));
        k13.end();
        
        k14 = y_values[j] - Data.M * (1-bin8_values[j]);
        sprintf(namme,"kkt14(%d)",int(j));
        kkt14.add(IloRange(env,-IloInfinity,k14,0,namme));
        k14.end();
        
        k15 = b_values[j] - Data.M * bin6_values[j];
        sprintf(namme,"kkt15(%d)",int(j));
        kkt15.add(IloRange(env,-IloInfinity,k15,0,namme));
        k15.end();
        
        k16 = 1 - y_values[j] - Data.M * (1-bin6_values[j]);
        sprintf(namme,"kkt16(%d)",int(j));
        kkt16.add(IloRange(env,-IloInfinity,k16,0,namme));
        k16.end();

    }

    model.add(kkt1);
    kkt1.end();
    model.add(kkt2);
    kkt2.end();
    model.add(kkt3);
    kkt3.end();
    model.add(kkt4); kkt4.end();
    model.add(kkt5); kkt5.end();
    model.add(kkt6); kkt6.end();
    model.add(kkt9); kkt9.end();
    model.add(kkt10); kkt10.end();
    model.add(kkt11); kkt11.end();
    model.add(kkt12); kkt12.end();
    model.add(kkt13); kkt13.end();
    model.add(kkt14); kkt14.end();
    model.add(kkt15); kkt15.end();
    model.add(kkt16); kkt16.end();

    ////////
    
    IloExpr kkt7(env);
    kkt7  = ubiom - Data.M * bin1 ;
    model.add(kkt7 <= 0);
    kkt7.end();
    
    IloExpr kkt8(env);
    kkt8  = flux_values[98] - Data.minbiomass - Data.M * (1-bin1) ;
    model.add(kkt8 <= 0);
    kkt8.end();
    
    IloExpr kkt17(env);
    kkt17  = c - Data.M * bin7 ;
    model.add(kkt17 <= 0);
    kkt17.end();
    
    IloExpr relax_kkt(env); // relaxation constraint
    IloExpr relax1_kkt(env);
    IloExpr kkt18(env);
    
    for (IloInt j = 0 ; j < nReac ; j++ ){
        relax_kkt += y_values[j];
        relax1_kkt += y_up_values[j];
    }
    kkt18 = (1+Data.ebs)*relax1_kkt - relax_kkt - Data.M * (1-bin7);
    model.add(kkt18 <= 0);
    relaxation1.end();
    relaxation.end();
    kkt18.end();
    
    ///////////////////////////////////////////    END OF CONSTRAINTS    ////////////////////////////////////////////
    

    IloCplex cplex(model);
    
    ////// to add starting point for binary knockout variables , z(j) //////////
    IloNumVarArray startVar(env);
    IloNumArray startVal(env);
    for (vector<int>::iterator i = Data.reaction_types.begin() ; i < Data.reaction_types.end() ; i++ ){
        if (*i == 2){
            startVar.add(z_values[i-Data.reaction_types.begin()]);
            startVal.add(1);
        }
    }
    cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto,  "MIPStart");
    startVal.end();
    startVar.end();
    
//    IloNumArray vals(env);
//    IloNumVarArray vars(env);
//    for (vector<int>::iterator i = Data.reaction_types.begin() ; i < Data.reaction_types.end() ; i++ ){
//        if (*i == 2){
//            vars.add(z_values[i-Data.reaction_types.begin()]);
//            vals.add(1);
//        }
//    }
//    cplex.setVectors(vals, 0, vars, 0, 0, 0);

    cplex.exportModel("/Users/meltemapaydin/Desktop/proom/moddel.lp");
        cplex.setParam(IloCplex::EpGap,0.01);
        cplex.solve();
        
        env.out() << endl << "Max succinate value: " << cplex.getObjValue() << endl;
    }
    catch(IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...){
        cerr << "Error: Unknown exception caught!!" << endl;
    }
    env.end();
    return 0;
}


