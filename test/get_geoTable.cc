#include <TFile.h>
#include <TMacro.h>
#include <TObjString.h>
#include <TList.h>
#include <iostream>
#include <fstream>

void test() {
    // Open ROOT file
    TFile *f = TFile::Open("/disk/data1/lze/ljuign/biasing/data_test_long/layer_10/B_Off_Concrete_Gammas_K40_0.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file.root" << std::endl;
        return;
    }

    // Get the TMacro
    TMacro *macro = (TMacro*)f->Get("geometryTable");
    if (!macro) {
        std::cerr << "geometryTable not found" << std::endl;
        return;
    }

    // Save to text
    std::ofstream out("geometryTable.txt");

    TList *lines = macro->GetListOfLines();
    if (lines) {
        for (int i = 0; i < lines->GetSize(); i++) {
            TObjString *obj = (TObjString*)lines->At(i);
            if (obj) out << obj->GetString().Data() << "\n";
        }
    }

    out.close();
    std::cout << "Saved geometryTable to geometryTable.txt" << std::endl;
}

def get_runbeamon_number(df: pd.DataFrame) -> int:
    """
    Given a pandas DataFrame from a runMacro TMacro,
    find the number after '/run/beamOn'.
    """
    # search for the row where first column is '/run/beamOn'
    mask = df.iloc[:,0] == '/run/beamOn'
    if mask.any():
        # grab the second column of that row (index 1)
        val = df.loc[mask, df.columns[1]].iloc[0]
        return int(val)
    else:
        raise ValueError("No '/run/beamOn' command found in DataFrame")
