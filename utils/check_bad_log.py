import os
import re

def search_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.stdout'):
                full_path = os.path.join(root, file)
                with open(full_path, 'r') as f:
                    content = f.read()
                    if re.search(r'Killed', content):
                        match = re.search(r'\w+\.root', content)
                        if match:
                            root_index = match.start()
                            before_root = content[:root_index].split()[-1]
                            print(f'{before_root} in {full_path}')

search_files('/eos/cms/store/group/phys_higgs/cmshgg/zhenxuan/custom_nanoAOD/log/UL17_TTJets/')
