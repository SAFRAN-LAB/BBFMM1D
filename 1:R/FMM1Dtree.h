#include<iostream>
#include<cmath>
#include<Eigen/Dense>
#include<vector>
const double PI	=	3.1415926535897932384;

// Class of LINE
class FMM1DLine{
    public:
    int lineNumber;  // Instead of box we can assume line for 1D
	int parentNumber;
	int childrenNumbers[2];
	int neighborNumbers[2];
	int innerNumbers[2]; // Inner interaction
	int outerNumbers[2]; // Outer interaction
	void printNodeDetails();
    //Constructor 
    FMM1DLine () {
	lineNumber		=	-1;
	parentNumber	=	-1;
	for (int l=0; l<2; ++l) {
		childrenNumbers[l]	=	-1;
        neighborNumbers[l]	=	-1;
        innerNumbers[l]		=	-1;
        outerNumbers[l]		=	-1;
		}
	}
        
	double center;   // Center of line
	Eigen::VectorXd multipoles;    
	Eigen::VectorXd locals;       
	//	The following will be stored only at the leaf nodes
	std::vector<double> chebNodes;
};

// Class of TREE
class FMM1DTree{
	public:
    int nLevels;			            		// Number of levels in the tree.
	int N;	                            		// Number of particles
    int nChebNodes;		                		// Number of Chebyshev nodes(still did not use)
    double L=1;                        			// Semi length of the simulation line
    int rank;				            		//	Rank of interaction, i.e., rank = nChebNodes*nChebNodes.
    std::vector<int> nLinesPerLevel;			//	Number of boxes per level in the tree
    std::vector<double> lineLength;				//	Line radius per level in the tree assuming the line at the root is [-1,1]
    std::vector<std::vector<FMM1DLine> > tree;	//	The tree storing all the information.
    void printTree();
    //	Chebyshev nodes
	std::vector<double> standardChebNodes1D;    
    std::vector<double> standardChebNodesChild;
	std::vector<double> leafChebNodes;
	// Different Operators
	Eigen::MatrixXd selfInteraction;		//	Needed only at the leaf level.
	Eigen::MatrixXd neighborInteraction[2];	//	Neighbor interaction only needed at the leaf level.
	Eigen::MatrixXd M2M[2];					//	Transfer from multipoles of 2 children to multipoles of parent.
	Eigen::MatrixXd L2L[2];					//	Transfer from locals of parent to locals of 2 children.
	Eigen::MatrixXd M2LInner[2];			//	M2L of inner interactions.
	Eigen::MatrixXd M2LOuter[2];			//	M2L of outer interactions.
	
	void printTreeDetails();
	void interactionDetails();
    
	// Constructor
	FMM1DTree(int nLevels, int nChebNodes) {
		this->nLevels			=	nLevels;
		this->nChebNodes		=	nChebNodes;
        this->rank				=	nChebNodes;
	    nLinesPerLevel.push_back(1);     //level 0
		lineLength.push_back(L);         //semi length
		for (int k=1; k<=nLevels; ++k) {
			nLinesPerLevel.push_back(2*nLinesPerLevel[k-1]);
			lineLength.push_back(0.5*lineLength[k-1]);
		}
		this->N					=	rank*nLinesPerLevel[nLevels];
	}
	
	// shifting Cheb nodes
	std::vector<double> shift_Cheb_Nodes(double xShift){
		std::vector<double> shiftedChebNodes;
			for (int k=0; k<rank; ++k) {
			double temp;
			temp	=	standardChebNodes1D[k]+2*xShift;
			shiftedChebNodes.push_back(temp);
		}
		return shiftedChebNodes;
	}
	// At leaf level
	std::vector<double> shift_Leaf_Cheb_Nodes(double xShift){
		std::vector<double> shiftedChebNodes;
			for (int k=0; k<rank; ++k) {
			double temp;
			temp	=	leafChebNodes[k]+xShift;
			shiftedChebNodes.push_back(temp);
		}
		return shiftedChebNodes;
	}
    
	// get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}
    // get_S as per paper
	double get_S(double x, double y, int n) {
		double S	=	0.0;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return (1.0/n)+2.0/n*S;
	}
    // set standard Cheb_Nodes
    void set_Standard_Cheb_Nodes() {
		for (int k=0; k<nChebNodes; ++k) {
			standardChebNodes1D.push_back(-cos((k+0.5)/nChebNodes*PI));
		}
		//	Left child, i.e., Child 0
        float temp;
		for (int j=0; j<rank; ++j) {
			temp	=	standardChebNodes1D[j];
			temp    =	0.5*temp-0.5;
			standardChebNodesChild.push_back(temp);
		}
		//	Right child, i.e., Child 1
		for (int j=0; j<rank; ++j) {
			temp	=	standardChebNodes1D[j];
			temp	=	0.5*temp+0.5;
			standardChebNodesChild.push_back(temp);
		}
		//std::cout<<"Here we can see all the ChebNodes for 2 children "<<std::endl;
		for(int j=0; j <standardChebNodesChild.size();j++){
			//std::cout << standardChebNodesChild[j]<<std::endl;
		}
    }
    void createTree() {
	//	First create root and add to tree
	FMM1DLine root;
	root.lineNumber		=	0;
	root.parentNumber	=	-1;
		
	for (int l=0; l<2; ++l) {
		root.childrenNumbers[l]	=	l;
        root.neighborNumbers[l]	=	-1;
        root.innerNumbers[l]	=	-1;
		root.outerNumbers[l]	=	-1;
	}
	
    std::vector<FMM1DLine> rootLevel;
	rootLevel.push_back(root);
	tree.push_back(rootLevel);

	for (int j=1; j<=nLevels; ++j) {
		std::vector<FMM1DLine> level;
		for (int k=0; k<nLinesPerLevel[j]; ++k) {
			FMM1DLine line;
			line.lineNumber		=	k;
			line.parentNumber	=	k/2;
			for (int l=0; l<2; ++l) {
				line.childrenNumbers[l]	=	2*k+l;
			}
			level.push_back(line);
		}
			tree.push_back(level);
		}
	}
	// Assigns the interactions for child 0 of a line
	void assign_Child0_Interaction(int j, int k){
		int nL	=	j+1;
		int nC	=	2*k;
		int nN;
		// Assigns siblings
		// ____**____|____N1____
		{	tree[nL][nC].neighborNumbers[1]	=	nC+1;
		}
		// Assign children of parent's 0th neighbor
		// ____I0____|____N0____
		nN	=	tree[j][k].neighborNumbers[0];
		if (nN != -1) {
			tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
			tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[1];		
		}
		// Assign children of parent's 1st neighbor
		// ____I1____|____O1____
		nN	=	tree[j][k].neighborNumbers[1];
		if (nN != -1) {
			tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
			tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[1];		
		}
		for(int l=0; l<2; l++){
			//std::cout << "Child 0 of line number ->" <<k << std::endl;
			//std::cout << "The " << l << "  neighbour is: " << tree[j][k].neighborNumbers[l] << std::endl;
			//std::cout << "The " << l << "  inner number is: " << tree[j][k].innerNumbers[l] << std::endl;
			//std::cout << "The " << l << "  outer number is: " << tree[j][k].outerNumbers[l] << std::endl;
		}
		//std::cout << "*******"<<std::endl;
	}
	// Assigns the interactions for child 1 of a line
	void assign_Child1_Interaction(int j, int k){
		int nL	=	j+1;
		int nC	=	2*k+1;
		int nN;
		// Assigns siblings
		// ____N0____|____**____
		{	tree[nL][nC].neighborNumbers[0]	=	nC-1;
		}
		// Assign children of parent's 0th neighbor
		// ____O0____|____I0____
		nN	=	tree[j][k].neighborNumbers[0];
		if (nN != -1) {
			tree[nL][nC].outerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
			tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[1];		
		}
		// Assign children of parent's 1st neighbor
		// ____N1____|____I1____
		nN	=	tree[j][k].neighborNumbers[1];
		if (nN != -1) {
			tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[0];
			tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[1];		
		}
		for(int l=0; l<2; l++){
			//std::cout << "Child 1 of line number ->" <<k << std::endl;
			//std::cout << "The " << l << "  neighbour is: " << tree[j][k].neighborNumbers[l] << std::endl;
			//std::cout << "The " << l << "  inner number is: " << tree[j][k].innerNumbers[l] << std::endl;
			//std::cout << "The " << l << "  outer number is: " << tree[j][k].outerNumbers[l] << std::endl;
		}
		//std::cout << "*******" << std::endl;
	}

	// Assigns the interactions for the children of a line
	void assign_Line_Interactions(int j, int k) {
		assign_Child0_Interaction(j,k);
		assign_Child1_Interaction(j,k);
	}

	// Assigns the interactions for the children all lines at a given level
	void assign_Level_Interactions(int j) {
		for (int k=0; k<nLinesPerLevel[j]; ++k) {	
			assign_Line_Interactions(j,k);
		}
	}

	// Assigns the interactions for the children all lines in the tree
	void assign_Tree_Interactions(){
		for (int j=0; j<nLevels; ++j){
			//std::cout << "At level " << j << std::endl;
			assign_Level_Interactions(j);
		}
	}
	void details_Of_Interactions(){
		for (int j=0; j<nLevels; j++){
			for (int k=0; k<nLinesPerLevel[j]; ++k){
				//std::cout << "interaction list of child 0 of line " << k <<" at level " << j <<std::endl;
				assign_Child0_Interaction(j,k);
				//std::cout << "interaction list of child 1 of line " << k <<" at level " << j <<std::endl;
				assign_Child1_Interaction(j,k);
				//std::cout << std::endl;
			}
			//std::cout << "~~~~~~~!! LEVEL END !!~~~~~~~" << std::endl;
		}
	}
	// Get different operators
	void get_Transfer_Matrix() {
		for (int l=0; l<2; ++l) {
			M2M[l]	=	Eigen::MatrixXd(rank,rank);
			for (int j=0; j<rank; ++j) {
				for (int k=0; k<rank; ++k) {
					M2M[l](j,k)	=	get_S(standardChebNodes1D[j], standardChebNodesChild[k+l*rank], nChebNodes);
				}
			}
			//std::cout << "Here is the M2M matrix for child "<< l << std::endl << std::endl << M2M[l] << std::endl;
			//std::cout << std::endl;
		}
		//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		for (int l=0; l<2; ++l) {
			L2L[l]	=	M2M[l].transpose();
			//std::cout << "Here is the L2L matrix for child "<< l << std::endl << std::endl << L2L[l] << std::endl;
			//std::cout << std::endl;
		}
	}
	// Assign charge at leaf level
	void assign_Leaf_Charges() {
		////std::cout << "Charge at leaf level :" << std::endl << std::endl;
		for (int k=0; k<nLinesPerLevel[nLevels]; ++k) {
			tree[nLevels][k].multipoles	= 0.5*(Eigen::VectorXd::Ones(rank)+Eigen::VectorXd::Random(rank));
			////std::cout << tree[nLevels][k].multipoles << std::endl<<std::endl;
		}
	}

	// Assign Kernel K(x,y)=1/|x-y|;
	double getInteraction(double x, double y){
		double s;
		s = fabs(x-y);
		if(s < 1e-10){
			return 0.0;       // If two chebs are equal then return 0
		}
		else{
			return 1.0/s;
		}
	}
	// We will use this to construct M2LInner and M2LOuter matrices
	void obtain_Desired_Operator(std::vector<double>& shiftedChebNodes, Eigen::MatrixXd& T) {
		T	=	Eigen::MatrixXd(rank,rank);
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				T(i,j)	=	getInteraction(standardChebNodes1D[i], shiftedChebNodes[j]);
			}
		}
	}

	void obtain_Desired_Leaf_Operator(std::vector<double>& shiftedChebNodes, Eigen::MatrixXd& T) {
		T	=	Eigen::MatrixXd(rank,rank);
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				T(i,j)	=	getInteraction(leafChebNodes[i], shiftedChebNodes[j]);
			}
		}

	}
	//	Assemble FMM Operators
	void assemble_Operators_FMM() {
		std::vector<double> shiftedChebNodes;
		//	Assemble Outer Interactions
		for (int l=0; l<2; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(6*l-3);
			//std::cout << "M2LOuter for " << l <<std::endl;
			obtain_Desired_Operator(shiftedChebNodes, M2LOuter[l]);
			//std::cout << M2LOuter[l] << std::endl;
			//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		}
		//	Assemble Inner Interactions
		for (int l=0; l<2; ++l) {
			shiftedChebNodes	=	shift_Cheb_Nodes(4*l-2);
			//std::cout << "M2LInner for " << l <<std::endl;
			obtain_Desired_Operator(shiftedChebNodes, M2LInner[l]);
			//std::cout << M2LInner[l] << std::endl;
			//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		}
		
		// Assigning Leaf Cheb
		for (int k=0; k<rank; ++k) {
			double temp;
			temp	=	lineLength[nLevels]*standardChebNodes1D[k];
			leafChebNodes.push_back(temp);
		}
		// Assemble Neighbor Interactions
		double neighborDistance	=	2.0*lineLength[nLevels];
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(-neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[0]);
			//std::cout << "N0 Matrix " << std::endl;
			//std::cout << neighborInteraction[0] << std::endl;
			//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

		}
		{
			shiftedChebNodes	=	shift_Leaf_Cheb_Nodes(neighborDistance);
			obtain_Desired_Leaf_Operator(shiftedChebNodes, neighborInteraction[1]);
			//std::cout << "N1 Matrix " << std::endl;
			//std::cout << neighborInteraction[1] << std::endl;
			//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		}
		// Assemble Self Interactions
		{
			obtain_Desired_Leaf_Operator(leafChebNodes, selfInteraction);
			//std::cout << "The self interaction Matrix is " << std::endl << std::endl;
			//std::cout << selfInteraction << std::endl;
			//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		}
	}
	// Assign Center location
	void assign_Center_Location() {
		int J, K;
		tree[0][0].center	=	0.0;
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*lineLength[j];
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				K	=	2*k;
				tree[J][K].center	=	tree[j][k].center-shift;
				tree[J][K+1].center	=	tree[j][k].center+shift;
			}
		}
		// Need to store Cheb nodes at leaf level
		for (int k=0; k<nLinesPerLevel[nLevels]; ++k) {
			tree[nLevels][k].chebNodes	=	shift_Leaf_Cheb_Nodes(tree[nLevels][k].center);
		}
	}

	// Evaluate M2M (upward pass)
	void evaluate_All_M2M() {
		for (int j=nLevels-1; j>1; --j) {
			int J	=	j+1;
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				int K	=	2*k;
				tree[j][k].multipoles	=	M2M[0]*tree[J][K].multipoles+M2M[1]*tree[J][K+1].multipoles;
			////std::cout << "Multipole at level  " << j << "  and line number " << k << std::endl << tree[j][k].multipoles << std::endl << std::endl;
			}
			////std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		}
	}
	// Evaluate M2L 
	void evaluate_All_M2L() {
		for (int j=2; j<=nLevels; ++j) {
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				tree[j][k].locals	=	Eigen::VectorXd::Zero(rank);
					// Inner well-separated clusters
					for (int l=0; l<2; ++l) {
						int nInner	=	tree[j][k].innerNumbers[l];
						//std::cout << "Inner Number " << std::endl;
						//std::cout << nInner << std::endl;
						if (nInner>-1) {
							tree[j][k].locals+=M2LInner[l]*tree[j][nInner].multipoles;
						}
					}
					// Outer well-separated clusters
					for (int l=0; l<2; ++l) {
						int nOuter	=	tree[j][k].outerNumbers[l];
						//std::cout << "Outer Number " << std::endl;
						//std::cout << nOuter << std::endl;
						if (nOuter>-1) {
							tree[j][k].locals+=M2LOuter[l]*tree[j][nOuter].multipoles;
						}
					}
				tree[j][k].locals/=lineLength[j];
			}
			//std::cout << "==========" << std::endl;
		}
	}
	// Evaluate L2L (Downward Pass)
	void evaluate_All_L2L() {
		for (int j=2; j<nLevels; ++j) {
			int J	=	j+1;
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				int K	=	2*k;
				tree[J][K].locals+=L2L[0]*tree[j][k].locals;
				tree[J][K+1].locals+=L2L[1]*tree[j][k].locals;
				// //std::cout << "At level " << j << " and line number " << k << " the local is " << std::endl << tree[j][k].locals << std::endl << std::endl;
			}
			// //std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
			// //std::cout << std::endl;
		}
	}
	// Evaluate nearby interaction at leaf level
	void evaluate_Leaf() {
		if (nLevels <2) {
			for (int k=0; k<nLinesPerLevel[nLevels]; ++k) {
				tree[nLevels][k].locals	=	Eigen::VectorXd::Zero(rank);
			}
		}
		for (int k=0; k<nLinesPerLevel[nLevels]; ++k) {
			for (int l=0; l<2; ++l) {
				int nNeighbor	=	tree[nLevels][k].neighborNumbers[l];
				if (nNeighbor > -1) {
					tree[nLevels][k].locals+=neighborInteraction[l]*tree[nLevels][nNeighbor].multipoles;
				}
			}
			tree[nLevels][k].locals+=selfInteraction*tree[nLevels][k].multipoles;
			//std::cout << "The required interaction using Blackbox FMM algo at line(nodes) " << k <<" is " << std::endl << tree[nLevels][k].locals << std::endl;
			//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		}
			//std::cout << std::endl;
			//std::cout << "Error checking ..." << std::endl << std::endl;
	}

	// Check error at a particular line
		void perform_Error_Check(int nLine) {
			Eigen::VectorXd potential	=	Eigen::VectorXd::Zero(rank);
			for (int l1=0; l1<rank; ++l1) {
				for (int k=0; k<nLinesPerLevel[nLevels]; ++k) {
					for (int l2=0; l2<rank; ++l2) {
						potential(l1)+=getInteraction(tree[nLevels][nLine].chebNodes[l1], tree[nLevels][k].chebNodes[l2])*(tree[nLevels][k].multipoles(l2));
					}
				}
			}
			//std::cout << "The EXACT potential at line(node) number "<< nLine << std::endl;
			//std::cout << potential << std::endl;  
			//std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
			Eigen::VectorXd error(rank);
			for (int k=0; k<rank; ++k) {
				error(k)	=	fabs((potential-tree[nLevels][nLine].locals)(k)/potential(k));
			}
			std::cout << "Error at line number " << nLine << " is" << std::endl;
			std::cout << error.maxCoeff();
		}
};


/*-------------------------------This section is only for debugging purpose-----------------------------------*/
// print node and tree details
void FMM1DLine::printNodeDetails(){	
    for(int i=0; i<2; i++){
		//std::cout << childrenNumbers[i] << std::endl;
	}
	 //std::cout << "~~~~~~~~~~~~~~~~~~~~~" <<std::endl;
}

// Prints details of all the lines in the tree:
void FMM1DTree::printTreeDetails(){
    for(int j = 0; j <= nLevels; j++){
		//std::cout << "At level " << j << std::endl<< std::endl;
        for(int k = 0; k < nLinesPerLevel[j]; k++){
			//std::cout << "Parent Number "<< std::endl;
			//std::cout << tree[j][k].parentNumber << std::endl;
			//std::cout << "line Number "<< std::endl;
			//std::cout << tree[j][k].lineNumber << std::endl;
			//std::cout << "Child Number "<< std::endl;
            tree[j][k].printNodeDetails();
        }
		//std::cout << "___!! LEVEL END !!___" << std::endl;
        //std::cout << std::endl << std::endl;
    }
}