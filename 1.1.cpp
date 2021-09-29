#include "ParseLatexFunc.h"
int main()
{
    string texStr(getLatexFileTemplate()),
           filename;
    ofstream log("log.txt");

    cout << "Enter fileName\n";
    cin >> filename;
    try{
        parseFile(filename, filename + "_tex.tex", filename + "_log.txt");
    }
    catch(exception& e){
        cout << e.what() << endl;
    }
    return 0;
}