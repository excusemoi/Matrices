#ifndef PARSELATEXFUNC_H
#define PARSELATEXFUNC_H
#include <iostream>
#include <fstream>
#include <math.h>
#include "SqMatrix.h"
#include <vector>
#include <stack>
using namespace std;

string getLatexFileTemplate()
{
    string latexFileTemplate(
           "\\documentclass[12pt]{article}\n"
           "\\usepackage{amsmath}\n"
           "\\title{Untitled Document}\n"
           "\\author{Valya Blin}\n"
           "\\date{\\today}\n"
           "\\begin{document}\n%\n"
           "\\end{document}\n"
           );
    return latexFileTemplate;
}

int isOperation(string operation)
{
    if(operation == "*")
        return 1;
    else if(operation == "+")
        return 1;
    else if(operation == "-")
        return 1;
    else if(operation == "==")
        return 1;
    else if(operation == "!=")
        return 1;
    return 0;
}

int isOperationCh(char ch)
{
    if(ch=='+'||ch=='-'||ch=='='||ch=='*'||ch=='!')
        return 1;
    return 0;
}

int isOperationStr(string str)
{
    if(str=="+"||str=="-"||str== "=="||str=="*"||str=="!=")
        return 1;
    return 0;
}

int isCharValid(char ch)
{
    if(isOperationCh(ch)||isalnum(ch)||ch=='.'||ch=='('||ch==')'||ch=='['||ch==']')
        return 1;
    return 0;
}

int isExpressionValid(string expression)
{
    stack<char> brackets;
    int length,
        exitCode = 0;
    for(int i = 0; i < expression.length(); i++){
        if(!isCharValid(expression[i]))
            throw invalid_argument("Incorrect characters in string");
        else if(expression[i] == '[' || expression[i] == '(')
            brackets.push(expression[i]);
        else if(expression[i] == ']')
            if(brackets.top() != '(')
                brackets.pop();
            else exitCode = -1;
        else if(expression[i] == ')')
            if(brackets.top() != '[')
                brackets.pop();
            else exitCode = -1;
        
        if(exitCode == -1)
            return 0;
    }
    return brackets.empty();
}

string eraseSpaces(string inputString)
{
    size_t index=0;
    while((index=inputString.find(' '))!=std::string::npos)
        inputString.erase(index, 1);
    return inputString;
}

string getDetTex(SqMatrix matr) 
{
    string matrix("\n$\\begin{vmatrix}\n");
    for(int i = 0; i < matr.getOrder(); i++)
        for(int j = 0; j < matr.getOrder(); j++){
            matrix += toString((matr)[i][j]);
            if(j != matr.getOrder() - 1)
                matrix += " & ";
            else
                if(i == matr.getOrder() - 1)
                    matrix += "\n\\end{vmatrix}$\n%";
                else 
                    matrix += "\\\\\n";
        }
    return matrix;
}

string getTexArithmOperation(SqMatrix firstMatrix, SqMatrix secondMatrix, string operation)
{
    SqMatrix resultMatrix(0, 0, false);
    double resultDouble;
    string resultString;
    if(!isOperationStr(operation))
        throw invalid_argument("Incorrect operation");
    else if((firstMatrix.getOrder() != secondMatrix.getOrder()) && operation != "==" && operation != "!=")
        throw invalid_argument("Incorrect orders of matrices");
    else if(operation == "*")
        resultMatrix = (firstMatrix * secondMatrix);
    else if(operation == "+")
        resultMatrix = firstMatrix + secondMatrix;
    else if(operation == "-")
        resultMatrix = firstMatrix - secondMatrix;
    else if(operation == "==")
        resultDouble = round(firstMatrix == secondMatrix);
    else if(operation == "!=")
        resultDouble = round(firstMatrix != secondMatrix);
    resultString += firstMatrix.convert();
    resultString += "\n$" + operation + "$";
    resultString += secondMatrix.convert();
    resultString += "\n$=$";
    resultString += (resultMatrix.getOrder() == 0)?("$"+toString(resultDouble)+"$"):(resultMatrix.convert());
    return resultString;
}

string getTexMultMatrixByNumber(SqMatrix matrix, double number)
{
    SqMatrix resultMatrix;
    string resultString;
    resultMatrix = matrix * number;
    resultString += matrix.convert();
    resultString += "\n$\\cdot$";
    resultString += "$" + toString(number) + "$";
    resultString += "$=$";
    resultString += resultMatrix.convert();
    return resultString;
}

string getTexDetReult(SqMatrix matrix)
{
    string resultString;
    double resultDouble;
    resultString += getDetTex(matrix);
    resultDouble = det(matrix);
    resultString += "\n$=$";
    resultString += toString(resultDouble);
    return resultString;
}

string getTexExpReult(SqMatrix matrix)
{
    string resultString;
    double resultDouble;
    SqMatrix resultMatrix = expMatrix(matrix, 4);
    resultString += "\nexp" + matrix.convert() +  "\n$=$\\mbox\n" + resultMatrix.convert() + "\n";
    return resultString;
}

string getTexTraceReult(SqMatrix matrix)
{
    string resultString;
    double resultDouble;
    resultDouble = traceMatrix(matrix);
    resultString += "\ntrace";
    resultString += matrix.convert() + "\n";
    resultString += "$=$\\mbox" + toString(resultDouble);
    return resultString;
}

string getTexTranspReult(SqMatrix matrix)
{
    string resultString;
    double resultDouble;
    SqMatrix resultMatrix;
    resultString += matrix.convert() + "\n$^{T}$" + "\n$=$\\mbox\n" + transp(resultMatrix = matrix).convert() + "\n";
    return resultString;
}

string getTexInverseReult(SqMatrix matrix)
{
    string resultString;
    double resultDouble;
    SqMatrix resultMatrix;
    if(fabs(det(matrix)) < EPS)
        throw invalid_argument("Inverse matrix of entered matrix doesn't exist");
    resultMatrix = inverseMatrix(matrix);
    resultString += matrix.convert() + "\n$^{-1}$" + "\n$=$\\mbox\n" + resultMatrix.convert() + "\n";
    return resultString;
}

string parseStringFromFile(string inputString)
{
    if(!isExpressionValid)
        throw invalid_argument("Incorrect brackets");
    inputString = eraseSpaces(inputString);
    int length = inputString.length(),
        nesting = 0,
        insideCommas = 0,
        outsideCommas = 0,
        order = 0;
    string buff = "",
           buff_ = "",
           operation = "",
           resultString = "";
    vector<int> rowsOrder;
    vector<double> matrixValues;
    vector<SqMatrix> matrices;
    SqMatrix resultMatrix(0, 0, false);
    double resultDouble;
    for(int i = 0; i < length; i++){
        if(inputString[i] == '[')
            nesting++;
        else if(inputString[i] == ']'){
            nesting--;
            if(nesting == 1){
                rowsOrder.push_back(insideCommas);
                matrixValues.push_back(stod(buff));
                buff = "";
                insideCommas = 0;
            }
            else if(nesting == 0){
                if(!(fabs((order = (int)sqrt(matrixValues.size()))*(int)sqrt(matrixValues.size())-matrixValues.size())<EPS))
                    throw invalid_argument("Matrix is non-square\n");
                for(auto& it : rowsOrder)
                    if(it != outsideCommas)
                        throw invalid_argument("Incorrect commas position");
                SqMatrix matrix(order, 0, false);
                for(int i = 0, k = 0; i < order; i++)
                    for(int j = 0; j < order; j++)
                        matrix[i][j] = matrixValues[k++];
                matrices.push_back(matrix);
                rowsOrder.clear();
                matrixValues.clear();
                outsideCommas = 0;
                continue;
            }
        }

        if(nesting > 2)
            throw invalid_argument("Incorrect dimension of matrix");
        else if(inputString[i] == ',')
            if(nesting == 1)
                outsideCommas++;
            else if(nesting == 2){
                insideCommas++;
                matrixValues.push_back(stod(buff));
                buff = "";
            }

        if(isdigit(inputString[i]) || inputString[i]=='.'\
         ||(inputString[i]=='-'&&!isdigit(inputString[i-1]))&&isdigit(inputString[i+1]))
            if(!nesting)
                buff_.push_back(inputString[i]);
            else
                buff.push_back(inputString[i]);
        else if(isalnum(inputString[i]) || isOperationCh(inputString[i]))
            operation.push_back(inputString[i]);
    }

    if(matrices.size() == 2)
        resultString = getTexArithmOperation(matrices[0], matrices[1], operation);
    else if(matrices.size() == 1 && buff_ != "" && operation != "==" && operation != "!=")
        resultString = getTexMultMatrixByNumber(matrices[0], stod(buff_));
    else if(matrices.size() == 1 && buff_ == ""){
        if(operation == "det")
            resultString = getTexDetReult(matrices[0]);
        else if(operation == "exp")
            resultString = getTexExpReult(matrices[0]);
        else if(operation == "trace")
            resultString = getTexTraceReult(matrices[0]);
        else if(operation == "inverse")
            resultString = getTexInverseReult(matrices[0]);
        else if(operation == "transp")
            resultString = getTexTranspReult(matrices[0]);
        else 
            throw invalid_argument("Incorrect input");
    }
    else
        throw invalid_argument("Incorrect input");
    return resultString;
}

int parseFile(string fileName, string texFileName, string logFileName)
{
    ifstream file(fileName);
    ofstream tex(texFileName), log(logFileName);
    string strFile, texStr(getLatexFileTemplate()), parsed;
    vector<string> info;
    int stringCounter = 1;
    if(!file.is_open() || !tex.is_open() || !log.is_open())
        throw invalid_argument("File can't be opened");
    while(getline(file, strFile)){
        try{
            if(strFile.empty()){
                stringCounter++;
                continue;
            }
            parsed = "\n" + parseStringFromFile(strFile) + "\n";
            info.push_back(parsed);
        }
        catch(exception& e){
            log << "Error in string #" << stringCounter << ":\n\t" << strFile <<"\n\t^\n\t" << e.what() << endl;
        }
        stringCounter++;
    }
    for(int i = 0; i < info.size(); i++)
        texStr.insert(texStr.find_last_of('%'), info[i]);
    tex << texStr;
    file.close();
    tex.close();
    log.close();
    system(("pdflatex " + texFileName).c_str());
    return 1;
}
#endif