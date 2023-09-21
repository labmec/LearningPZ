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
//Cria a classe EMatid, que é do tipo enum, isso significa que ENone, EDomain e EBC são variáveis que não podem ter seu valor alterado diretamente, preservando-as em seu estado original relacionado a classe enum.
const int global_nthread = 16;
// Determina o valor máximo de threads para o matskl.SetNumThreads, isto é, 16.

TPZGeoMesh* CreateGMesh(int ndiv);
// Define a função CreateGMesh, que recebe um valor inteiro ndiv que é retornado como um pointer a TPZGeoMesh.
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
// Analogamente, // Define a função ReadMeshFromGmsh, que recebe uma string relacionada ao nome de um arquivo, que é retornado como um pointer a TPZGeoMesh.
void CreateBCs(TPZGeoMesh* gmesh);
// Cria a função CreateBCs que não retorna valor, cuja entrada é o pointer para TPZGeoMesh.
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas);
// CreateH1CMesh é uma função que retorna um pointer a TPZCompMesh, que tem os seguintes parâmetros:
//TPZGeoMesh* gmesh é um pointer para TPZGeoMesh
// const int pord, é uma constante da variável inteira pord e, portanto, não pode ser alterada dentro da função
//TElasticity3DAnalytic *elas é um pointer para TElasticity3DAnalytic
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
// Define a função do tipo void (Não retorna valor) SolveProblemDirect, que contém os seguintes parâmetros:
//TPZLinearAnalysis &an, é uma referência a função TPZLinearAnalysis denominada "an". Desse modo, quando essa referência, denominada "an", é alterada, o objeto original de onde ela foi retirada também será alterado, e não uma cópia desse objeto.
//TPZCompMesh *cmesh é um pointer para um objeto do tipo TPZCompMesh.
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
//É definida a função PrintResults do tipo void, cujos parâmetros são os mesmos já explicados na linha anterior.


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
// Associa uma forma de armazenar o local de memória de TElasticity3DAnalytic denominado "elas", em um novo objeto, devido ao operador new
    elas->fE = 250.;//206.8150271873455;
//Define o valor da variável fE, pertencente ao objeto cujo pointer é elas (que, pelo código, aparenta ser TElasticity3DAnaçytic), para 250 e, provavelmente, antes era para 206.8150271873455 e tal valor foi alterado para um comentário. 
    elas->fPoisson = 0.;
//Define o valor de fPoisson, pertencente a TElasticity3DAnalytic devido ao pointer "elas", para 0.
    elas->fProblemType = TElasticity3DAnalytic::EStretchx;
//Essa linha é relacionada ao valor de fProblemType, pertencente a TElasticity3DAnalytic devido ao pointer "elas". Desse modo, ele define um valor do tipo enum para a variável trabalhada.
    TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh,pord,elas);
// Define o pointer chamado cmeshH1 relacionado a classe TPZCompMesh, associando-o a função CreateH1CMesh, cujos parâmetros (já antes trabalhados) são: gmesh ,pord,elas.
    
    TPZLinearAnalysis an(cmeshH1);
// Nesse caso, a adição do termo "an" significa que estamos criando um objeto da classe "TPZLinearAnalysis", cujo nome é "an", e tal sintaxe é a de um constructor, inicializando com o objeto "an" com o argumento cmeshH1.
    SolveProblemDirect(an,cmeshH1);
// A função SolveProblemDirect recebe dois parâmetros - an e cmeshH1. an é o objeto de referência recém criado, relacionado a TPZLinearAnalysis, enquanto que cmeshH1 é um pointer a um objeto de TPZCompMesh.
    
    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
// Printa na tela "--------- PostProcess ---------"
    PrintResults(an,cmeshH1);
// Define a função PrintResults, que recebe os mesmos parâmetros recém explicados,
    
    // deleting stuff
    delete cmeshH1;
    delete gmesh;
// Como o comentário sugere, ambas as linhas acima são aplicadas para deletar os valores de cmeshH1 e gmesh, a fim de liberar espaço de memória de armazenamento, se livrando de dados que não serão mais utilizados.
        
    std::cout << "--------- Simulation finished ---------" << std::endl;
// Printa na tela "--------- Simulation finished ---------"
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
// Cria a função CreateH1CMesh que retornará um pointer associado a classe TPZCompMesh, que aceita os seguintes parâmetros: TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas.

    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
//
    const int dim = gmesh->Dimension();
//
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
// Essas sintaxes tem a função de "chamar" o pointer cmesh de modo a acessar membros pertencentes a ele e associar a diferentes funções, explicitadas acima, cada uma contendo um argumento refente a seu próprio código.
    
    // Domain elas mat
    const STATE E = elas->fE, nu = elas->fPoisson;
//
    TPZManVector<STATE> force = {0,0,0};
//
    TPZElasticity3D *mat = new TPZElasticity3D(EDomain, E, nu, force, 0., 0., 0.);
//
    mat->SetExactSol(elas->ExactSolution(), 2);
//
    mat->SetForcingFunction(elas->ForceFunc(), 4);
//
    cmesh->InsertMaterialObject(mat);
//
    
    // BC null mat
    TPZFMatrix<STATE> val1(3,3,0.);
//
    TPZManVector<STATE> val2(3,0.);
//
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
// Essa linha só define valores arbitrários para diversas variáveis, a fim de servirem como parâmetros a serem utilizados em outras funções.
    auto* BCCond0 = mat->CreateBC(mat, EBC, diri, val1, val2);
//
    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 4);
//
    cmesh->InsertMaterialObject(BCCond0);
// Essa sintaxe chama o pointer cmesh e o associa a uma função denominada InsertMaterialObject, cujo parâmetro de entrada é BCCond0
    
    // Constructs mesh
    cmesh->AutoBuild();
// Essa linha de código "chama' a função AutoBuild dentro do objeto ao qual cmesh é pointer"
    
    return cmesh;
// Define como valor a ser retornado a variável cmesh, após concluir a função.
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
// Declara a função SolveProblemDirect do tipo void (logo não retornará valor), com os seguintes parâmetros:
//TPZLinearAnalysis &an, cuja sintaxe demonstra ser uma referência ao TPZLinearAnalysis
//TPZCompMesh *cmesh, cuja sintaxe é
{

    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    // TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
//
    an.SetStructuralMatrix(matskl);
//
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
//
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
//
    an.SetSolver(step);
//
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
// printa a linha "--------- Assemble ---------"
    TPZSimpleTimer time_ass;
//
    an.Assemble();
//
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;
// printa para o usuário a linha "Total time = (Valor resultante do cálculo de time_ass.ReturnTimeDouble()/1000.)s"

    
    ///solves the system
    std::cout << "--------- Solve ---------" << std::endl;
// printa a linha "--------- Solve ---------"
    TPZSimpleTimer time_sol;
//
    an.Solve();
//
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
// printa para o usuário a linha "Total time = (Valor resultante do cálculo de time_sol.ReturnTimeDouble()/1000.)s"
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
//Define a função do tipo void chamada PrintResults, que recebe como parâmetros TPZLinearAnalysis &an e  TPZCompMesh *cmesh
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
// printa para o usuário a linha "--------- Post Process ---------"
    TPZSimpleTimer postProc("Post processing time");
// Declara uma variável chamada postProc, do tipo TPZSimpleTimer, chamando um construtor com uma string como argumento, igual a "Post processing time".
    const std::string plotfile = "postprocess";
// define a string plotfile, do tipo const, equivalente a uma string de texto "postprocess"
    constexpr int vtkRes{0};
// define a variável do tipo inteiro denominada vtkRes, do tipo constexpr, que significa que é uma expressão constante, ou seja,  vtkRes é um valor constante e não pode ser alterado. Ainda, {0} indica o valor associado a essa constante, e portanto não será alterado, com valor determinado na hora de compilação.
    
    TPZVec<std::string> fields = {
        // "ExactDisplacement",
        // "ExactStress",
        "Displacement",
        "Stress",
        "Strain",
    };
// Nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas ao deslocamento, deformação e tensão. Além disso, mais duas strings estão no código porém como estão comentadas, devem ter sido o nome antigo dado a elas.
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
// Essa linha de código declara uma variável chamada vtk do tipo auto, o que significa que o compilador irá deduzir o tipo que ela terá a depender do que ela é igual. No caso, ela é igual a função TPZVTKGenerator, de parâmetros cmesh, fields, plotfile, vtkRes.
    vtk.SetNThreads(global_nthread);
// Chama uma função "SetNThreads", de parâmetro "global_nthread" no objeto vtk. Desse modo, a função vtk realiza uma ação a partir da entrada fornecida que é a função SetNThreads, utilizando a estrutura de função do objeto vtk.
    vtk.Do();
// Chama uma função chamada "Do" no objeto "vtk". Basicamente, representa uma ação que está relacionada a função vtk.
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
// printa para o usuário a linha "Total time = (Valor resultante do cálculo de postProc.ReturnTimeDouble()/1000.)s"
    
    return;
}