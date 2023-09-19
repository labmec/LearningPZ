#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZNullMaterial.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <TPZSimpleTimer.h>
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include <Elasticity/TPZElasticity3D.h>
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "tpzchangeel.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

enum EMatid {ENone,EDomain,EBC};
const int global_nthread = 16;


TPZGeoMesh* CreateGMesh(int ndiv);
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);


int main() {
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    //printa o texto "--------- Starting simulation ---------", indicando que a simulação está se iniciando.
    const int pord = 1;
    //é declarada a constante pord, cujo valor inteiro atribuído é 1.
    int ndiv = 2;
    //é declarada a variável ndiv (representa o número de divisões de cada uma das dimensões da malha), cujo valor inteiro atribuído é 2.
    TPZGeoMesh* gmesh = CreateGMesh(ndiv);
    //criação da váriável gmesh (geometric mesh) cujo valor será um ponteiro de um objeto da classe TPZGeoMesh, atribuído através da função CreateGMesh(ndiv), que recebe o valor de ndiv como entrada. 
    //o que ocorre na função está detalhado mais a frente no código.    
    std::ofstream out("gmesh.vtk");
    //a variável out está sendo declarada como pertencente à classe std::ofstream (usada para lidar com arquivos de saída). o nome do arquivo que será aberto para escrita é "gmesh.vtk". 
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    //está sendo chamada a função PrintGMeshVTK da classe TPZVTKGeoMesh, cujos argumentos são "gmesh" e "out".
    //essa função formata os dados da malha geométrica "gmesh" para o formato VTK e os escreve no objeto "out", que representa o arquivo "gmesh.vtk"

    TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
    elas->fE = 250.;//206.8150271873455;
    elas->fPoisson = 0.;
    elas->fProblemType = TElasticity3DAnalytic::EStretchx;
    TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh,pord,elas);
    
    TPZLinearAnalysis an(cmeshH1);
    SolveProblemDirect(an,cmeshH1);
    
    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an,cmeshH1);
    
    // deleting stuff
    delete cmeshH1;
    delete gmesh;
        
    std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZGeoMesh* CreateGMesh(int ndiv) {
    //está sendo definida a função CreateGMesh que retornará um ponteiro para um objeto da classe TPZGeoMesh, a função recebe o argumento inteiro "ndiv".
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    //gmesh é o endereço de um novo objeto criado pertencente à classe TPZGeoMesh.
    
    MMeshType meshType = MMeshType::EHexahedral;
    //a variável meshType está recebendo o valor EHexahedral, indicando que os elementos da malha serão hexaédricos (cubos - 6 faces).
    int dim = 3;
    //variável "dim" é definida com valor inteiro 3 indicando que a malha será tridimensional.
    TPZManVector<REAL,3> minX = {-1,-1,-1};
    TPZManVector<REAL,3> maxX = {1,1,1};
    //Ambos os vetores "minX" e "maxX" possuem 3 elementos que indicam as coordenadas máximas e mínimas, respectivamente, da malha geométrica em cada uma das 3 dimensões.
    //Nesse caso, a caixa da malha geométrica vai de -1 a 1 nas 3 dimensões.
    //TPZManVector é uma classe utilizada para criar vetores desse tipo.
    int nMats = 2*dim+1;
    // número de materiais talvez.
    
    constexpr bool createBoundEls{true};
    //é atribuído "true" à váriavel booleana "createBoundEls", indicando provavelmente que condições de contorno serão inseridas na malha.
    TPZVec<int> matIds(nMats,EBC);
    //é criado o vetor inteiro matIds com "nMats" elemento(s) cujos valores serão todos "EBC".
    matIds[0] = EDomain;
    //o primeiro elemento (índice 0) do vetor matIds está sendo definido como Edomain.
    
    TPZManVector<int,3> ndivvec = {ndiv,ndiv,ndiv};
    //é criado o vetor inteiro ndivvec com 3 elementos, todos com valor ndiv, indicando que cada uma das 3 dimensões da malha será dividida em ndiv (nesse caso 2) partes, ou seja, 2 elementos.
    gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,matIds, ndivvec, meshType,createBoundEls);
    //é chamada a função CreateGeoMeshOnGrid para criar finalmente a malha geométrica com os parâmetros específicados anteriormente. o resultado é atribuído à variável gmesh.

    return gmesh;
    //a função retorna o ponteiro gmesh que aponta para a malha criada.
}


TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas) {
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    // Domain elas mat
    const STATE E = elas->fE, nu = elas->fPoisson;
    TPZManVector<STATE> force = {0,0,0};
    TPZElasticity3D *mat = new TPZElasticity3D(EDomain, E, nu, force, 0., 0., 0.);
    mat->SetExactSol(elas->ExactSolution(), 2);
    mat->SetForcingFunction(elas->ForceFunc(), 4);
    cmesh->InsertMaterialObject(mat);
    
    // BC null mat
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 4);
    cmesh->InsertMaterialObject(BCCond0);
    
    // Constructs mesh
    cmesh->AutoBuild();
    
    return cmesh;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{

    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    // TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    ///solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess";
    constexpr int vtkRes{0};
    
    TPZVec<std::string> fields = {
        // "ExactDisplacement",
        // "ExactStress",
        "Displacement",
        "Stress",
        "Strain",
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}