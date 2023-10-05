#include "TFile.h"
#include "TTree.h"

void func_add_weight_central(const char* inputFilename, const char* outputFilename) {
    // 打开原始ROOT文件
    TFile* inputFile = new TFile(inputFilename, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file: " << inputFilename << std::endl;
        return;
    }

    // 获取原始TTree
    TTree* inputTree = (TTree*)inputFile->Get("cat1");
    if (!inputTree) {
        std::cerr << "Error reading input tree." << std::endl;
        return;
    }

    // 创建新的ROOT文件
    TFile* outputFile = new TFile(outputFilename, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFilename << std::endl;
        return;
    }

    // 创建新的TTree，并复制所有原始分支
    TTree* outputTree = inputTree->CloneTree(0);

    // 获取原始分支
    float weight, weight_central;


    // 创建新的分支
    TBranch* minBranch = outputTree->Branch("weight_central", &weight_central, "weight_central/F");

    // 循环所有事件，计算最小值和最大值，并填充到新的分支
    for (int i = 0; i < inputTree->GetEntries(); i++) {
        inputTree->GetEntry(i);
        weight_central = 1;
        // minBranch->Fill();
        // maxBranch->Fill();
        outputTree->Fill();
    }

    // 写入新的ROOT文件
    outputFile->Write();
    outputFile->Close();
    inputFile->Close();

    // 释放内存
    delete outputFile;
    delete inputFile;
}
void add_weight_central() {
    // const char* inputFilename = "/eos/user/z/zhenxuan/hhwwgg_root/data_UL17cat1.root";
    // const char* outputFilename = "/eos/user/z/zhenxuan/hhwwgg_root/data_UL17_cat1.root";
    // addMinAndMaxPhotonID(inputFilename, outputFilename);
    const char* inputFilename = "/eos/user/z/zhenxuan/hhwwgg_root/combined_data_cat1_minmax.root";
    const char* outputFilename = "/eos/user/z/zhenxuan/hhwwgg_root/combined_data_cat1_minmax_weight_central.root";
    func_add_weight_central(inputFilename, outputFilename);

    return 0;
}
