#include <iostream>
#include <chrono>
#include <ctime>
#include "FMM1Dlogtree.hpp"
using namespace std;
int main(int argc, char* argv[]){
	int nLevels		=	atoi(argv[1]);
	int nChebNodes	=	atoi(argv[2]);
    int line        =   atoi(argv[3])            ;
	int L	        =    1;	
    std::chrono::time_point<std::chrono::system_clock> start, end;
    FMM1DTree t1(nLevels,nChebNodes);
    start = std::chrono::system_clock::now();
    t1.set_Standard_Cheb_Nodes();
    t1.createTree();
    t1.printTreeDetails();
    // cout << "### Interaction details ###" <<endl<<endl;
    //t1.details_Of_Interactions();
    t1.assign_Tree_Interactions();
    t1.get_Transfer_Matrix();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_tree = end - start;

    start = std::chrono::system_clock::now();
    t1.assemble_Operators_FMM();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_assemble = end - start;
    
    start = std::chrono::system_clock::now();
    t1.assign_Leaf_Charges();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_charge = end - start;
    
    start = std::chrono::system_clock::now();
    t1.evaluate_All_M2M();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_M2M = end - start;
   
    start = std::chrono::system_clock::now();
    t1.evaluate_All_M2L();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_M2L = end - start;
    
    start = std::chrono::system_clock::now();
    t1.evaluate_All_L2L();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_L2L = end - start;
    
    start = std::chrono::system_clock::now();
    t1.evaluate_Leaf();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_leaf = end - start;
    cout <<"Performing ERROR checking..." << endl << endl;
    start = std::chrono::system_clock::now();
    t1.assign_Center_Location();
    t1.perform_Error_Check(line);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_error = end - start;
    cout << endl <<endl;

    cout <<"********** Summary **********" << endl << endl;
    cout << "The number of particles taken: " << pow(2,nLevels)*nChebNodes << endl << endl;
    cout << "The number of Cheb Nodes taken: " << nChebNodes << endl << endl;
    cout << "You wanna check error at line: " << line << endl << endl;
    cout << "elapsed time to create tree: " << elapsed_tree.count() << "s\n" << endl;
    cout << "elapsed time for assigning charges at leaf level: " << elapsed_charge.count() << "s\n" << endl;
    cout << "elapsed time for assembling operators: " << elapsed_assemble.count() << "s\n" << endl;
    cout << "elapsed time for M2M: " << elapsed_M2M.count() << "s\n" << endl;
    cout << "elapsed time for M2L: " << elapsed_M2L.count() << "s\n" << endl;
    cout << "elapsed time for L2L: " << elapsed_L2L.count() << "s\n" << endl;
    cout << "elapsed time for Nearby interaction at leaf level: " << elapsed_leaf.count() << "s\n" << endl;
    cout << "*** Total time taken nearly -> " << elapsed_L2L.count()+elapsed_M2L.count()+elapsed_M2M.count()+elapsed_assemble.count()+elapsed_charge.count()+elapsed_tree.count()+elapsed_leaf.count() << "s\n" << endl;
    cout << "elapsed time for calculating error at a given line: " << elapsed_leaf.count() << "s\n" << endl;
    
}
