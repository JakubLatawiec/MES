#pragma once

#include <array>
#include <vector>

#include "Jacobian.h"
#include "ElementUniv.h"
#include "Node.h"
#include "../utils/Gauss.h"
#include "Surface.h"

class Element
{
private:
	static int m_IPC;
	static int m_IPC2D;
	static int m_SFC;
	static std::unique_ptr<ElementUniv> m_ElementUniv;
	static const ElementUniv& getElementUniv();

private:
	std::array<int, 4> m_NodesID{};
	std::vector<Jacobian> m_JacobianMatrixes{};
	std::vector<Matrix> m_StiffnessMatrixes{};
	Matrix m_StiffnessMatrix{};
	std::array<Matrix, 4> m_HbcMatrixes{};
	Matrix m_PVector{};
	std::vector<Matrix> m_CMatrixes{};
	Matrix m_CMatrix{};


public:
	//Constructors
	Element() = default;
	
	//Calculation methods
	void CalcJacobians(const std::vector<Node>& nodes);
	void CalcStiffnessMatrixes(double conductivity);
	void CalcHbcMatrixes(const std::vector<Node>& nodes, double alpha);
	void CalcPVector(const std::vector<Node>& nodes, double alpha, double Tot);
	void CalcCMatrix(double c, double rho);

	static void Initialize(int ipc, int sfc);

	//Debug methods
	static void PrintElementUniv();
	void PrintJacobianMatrixes();
	void PrintStiffnessMatrixes();
	void PrintStiffnessMatrix();
	void PrintNodesID();
	void PrintHbcMatrixes();
	void PrintPVector();
	void PrintCMatrixes();
	void PrintCMatrix();

	//Setters
	void setNodesID(size_t index, int id);

	//Getters
	const std::array<int, 4>& getNodesID() const;
	const Matrix& getStifnessMatrix() const;
	const Matrix& getPVector() const;
	const Matrix& getCMatrix() const;
};