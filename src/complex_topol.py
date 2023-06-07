import sys
import argparse
import textwrap


#Função para encontrar numero da linha que contém o argumemnto desejado
def get_line(txt):
    file = open('topol.top', encoding='utf8')
    for line_num, value in enumerate(file, 1): #le todas as linhas e valores do file começando como 1
        if txt in value:                       #qdo encontra o 'txt' retorna o valor da linha
            return line_num


if len(sys.argv) == 2 and sys.argv[1] == '-h':


    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                          Este script deve ser rodado no diretorio em analise
                                          Ele cria o arquivo complex.gro e faz as adequações ao topol.top
                                          Usado para proteina e ligante ou proteina ligante e cofator
                                          sintaxe -> python3.9 complex.py prot_name lig_name cof_name
                                          as extensoes DEVEM ser adicionadas ao comando
                                          '''))

    args = parser.parse_args()


if len(sys.argv) > 4:
    print('ERROR: Please put protein ligand and cofactor names or -h for help. ')


#-----------------------------------------------------------------------------------------------------

elif len(sys.argv) == 3:
    prot = sys.argv[1]
    lig = sys.argv[2]
    prot2 = prot.rsplit('_', 1)[0]
    lig2 = lig.rsplit('.', 1)[0]

    # Abre arquivo e lê prot
    file = open(str(prot), 'r')
    prot = file.readlines()

    # Criar complex, parte1 -> prot sem box vector (only_prot)
    file = open('complex.gro', 'a+')
    only_prot = prot[2:-1]
    file.writelines([item for item in only_prot])
    box_vector = prot[-1]
    file.close()

    # incluir ligante sem cabeçalho nem box vector -> incluir box vector
    # Abre arquivo e lê lig e complex
    file = open(str(lig), 'r')
    lig = file.readlines()
    lig_only = lig[2:-1]

    file = open('complex.gro', 'a+')
    file.writelines([item for item in lig_only])
    file.writelines([item for item in box_vector])
    file.close()

    # Verificar o numero de linhas (cont) e a linha2 - como esta sem cabeçalho o numero de moleculas é cont -1
    file = open('complex.gro', 'r')
    line = file.readlines()
    cont = len(line)
    file.close()

    # Substituir o numero de linhas e o nome inicial
    file = open('complex.gro', 'w')
    file.write('complex' + '\n')
    file.write(str(cont - 1) + '\n')
    file = open('complex.gro', 'a')
    file.writelines(line)
    file.close()

#--------------------------------------------------------------------------------------------

    x = str('#include "amber99sb.ff/forcefield.itp"')
    y = str('; Include Position restraint file')
    z = str('#mols')

    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(x) + 1, f'; Include ligand parameters\n#include "./{lig2}.prm"\n\n')
    file.seek(0)
    file.writelines(line)

    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(x) + 4, f'; Include protein topology\n#include "./protein.itp"\n\n')
    file.seek(0)
    file.writelines(line)

    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(y) + 3, f'\n; Include ligand topology\n#include "./{lig2}.itp"\n')
    file.seek(0)
    file.writelines(line)

    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(z) + 1, f'{prot2.capitalize()}             1\n{lig2.upper()}             1\n')
    file.seek(0)
    file.writelines(line)

#---------------------------------------------------------------------------------------------

elif len(sys.argv) == 4:
    prot = sys.argv[1]
    lig = sys.argv[2]
    cof = sys.argv[3]
    prot2 = prot.rsplit('_', 1)[0]
    lig2 = lig.rsplit('.', 1)[0]
    cof2 = cof.rsplit('.', 1)[0]

    # Abre arquivo e lê prot
    file = open(str(prot), 'r')
    prot = file.readlines()

    #Criar complex, parte1 -> prot sem box vector (only_prot)
    file = open('complex.gro', 'a+')
    only_prot = prot[2:-1]
    file.writelines([item for item in only_prot])
    box_vector = prot[-1]
    file.close()

    #incluir ligante sem cabeçalho nem box vector -> incluir box vector
    #Abre arquivo e lê lig e complex
    file = open(str(lig), 'r')
    lig = file.readlines()
    lig_only = lig[2:-1]

    file = open(str(cof), 'r')
    cof = file.readlines()
    cof_only = cof[2:-1]

    file = open('complex.gro', 'a+')
    file.writelines([item for item in lig_only])
    file.writelines([item for item in cof_only])
    file.writelines([item for item in box_vector])
    file.close()

    #Verificar o numero de linhas (cont) e a linha2 - como esta sem cabeçalho o numero de moleculas é cont -1
    file = open('complex.gro', 'r')
    line = file.readlines()
    cont = len(line)
    file.close()

    #Substituir o numero de linhas e o nome inicial
    file = open('complex.gro', 'w')
    file.write('complex' + '\n')
    file.write(str(cont -1)+ '\n')
    file = open('complex.gro', 'a')
    file.writelines(line)
    file.close()

#-----------------------------------------------------------------


    x = str('#include "amber99sb.ff/forcefield.itp')
    y = str('; Include Position restraint file')
    z = str('#mols')

    # Abrir o topol.top -> Ler as linhas a add o argumento .prm
    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(x), f'#include "./protein.itp"\n')
    file.seek(0)
    file.writelines(line)

    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(x)+1, f'#include "./{lig2}.rtp"\n#include "./{cof2}.rtp"\n')
    file.seek(0)
    file.writelines(line)

    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(y) + 3, f'\n#ifdef POSRES_LIG\n#include "./posre_{lig2}.itp"\n#endif\n')
    file.seek(0)
    file.writelines(line)

    # Abrir o topol.top -> Ler as linhas a add o argumento .itp
    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(y) + 10, f'\n; Include ligand topology\n#include "./{lig2}.itp"\n#include "./{cof2}.itp"\n')
    file.seek(0)
    file.writelines(line)

    # Abrir o topol.top -> Ler as linhas a add o argumento numero de moleculas finais
    file = open('topol.top', 'r+', encoding='utf8')
    line = file.readlines()
    line.insert(get_line(z) + 1, f'{prot2.capitalize()}             1\n{lig2.upper()}             1\n{cof2.upper()}             1\n')
    file.seek(0)
    file.writelines(line)
