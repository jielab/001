import re
import os
import sys
import pandas as pd


def extract_table(read_text: str, output_dir: str = './res'):
    type_list = '新增确诊病例|境外输入病例|本土病例|疑似病例|治愈出院病例|现有治愈出院病' \
                '解除医学观察的密切接触者|重症病例|确诊病例|' \
                '累计治愈出院病|死亡病例|密切接触者|医学观察的密切接触者|无症状感染者|当日转为确诊病例|' \
                '当日解除医学观察|尚在医学观察的无症状感染者|由无症状感染者转为确诊病例|出院|死亡'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def get_postion(positions: list, type_collect: str):
        if not positions:
            position = '31个省（自治区、直辖市）和新疆生产建设兵团' if type_collect != '现有（港澳台）' else '港澳台地区'
        else:
            position = positions[i - 1]
        return position

    def split_sentence(paragraph: str):
        semicolon_pos = [i.start() for i in re.finditer('；', paragraph)]
        bracket_left_pos = [i.start() for i in re.finditer('（', paragraph)]
        bracket_right_pos = [i.start() for i in re.finditer('）', paragraph)]
        assert len(bracket_left_pos) == len(bracket_right_pos)
        in_bracket = [False for _ in range(len(semicolon_pos))]
        for idx, pos in enumerate(semicolon_pos):
            for i in range(len(bracket_left_pos)):
                if bracket_left_pos[i] < pos < bracket_right_pos[i]:
                    in_bracket[idx] = True
        split_pos = [pos for idx, pos in enumerate(semicolon_pos) if not in_bracket[idx]]
        split_pos = [-1] + split_pos + [len(paragraph)]
        return [paragraph[split_pos[idx] + 1:split_pos[idx + 1]] for idx in range(len(split_pos) - 1)]

    with open(read_text, 'r', encoding='UTF-8') as f:
        text = f.readlines()

    while '\n' in text:
        text.remove('\n')

    result = pd.DataFrame()
    date = re.search(r'\d{4}-\d{2}-\d{2}', ''.join(text))[0].replace('-', '_')

    for row in text:
        type_collect = '现有' if '现有' in row else '新增' if '新增' in row else '现有（港澳台）' if '港澳台' in row else '其它'
        row = re.sub('。|，含', '；', row)
        sentences = split_sentence(row)
        abroad = True if row.startswith('境外输入') else False
        for sentence in sentences:
            phrases = re.split('[，、；（）]', sentence)

            types = []
            positions = []
            numbers = []
            i = 0
            for phrase in phrases:
                if re.findall(r'\d+[例人]', phrase):
                    prefix = re.findall(r'[\D]*(?=\d+[例人])', phrase)[0]
                    number = int(re.findall(r'\d+(?=[例人])', phrase)[0])
                    if re.findall(type_list, phrase):
                        type_ = re.findall(type_list, phrase)[0]
                        position = get_postion(positions, type_collect)
                        type_ = '境外输入' + type_ if abroad else type_
                    else:
                        position = prefix
                        if '其中' in position:
                            position = re.sub('其中', '', position)

                        if types:
                            if type_collect == '现有（港澳台）':
                                i = 1
                                while types[i - i] == '出院' or types[i - i] == '死亡':
                                    i += 1
                                type_ = types[i - i]
                            else:
                                type_ = types[i - 1]
                        else:
                            if '其中' in sentence:
                                type_ = result.iloc[-1, 1]
                            else:
                                raise Exception('Error')
                elif re.findall('均在', phrase):
                    suffix = re.findall(r'(?<=均在)[\D]+', phrase)
                    number = numbers[i - 1]
                    position = suffix[0]
                    type_ = types[i - 1]
                elif re.findall('均为', phrase):
                    suffix = re.findall(r'(?<=均为)[\D]+', phrase)
                    position = get_postion(positions, type_collect)
                    type_ = suffix[0]
                    type_ = '境外输入' + type_ if abroad else type_
                    number = numbers[i - 1]
                elif re.findall('无', phrase):
                    suffix = re.findall(r'(?<=无)[\D]+', phrase)
                    position = get_postion(positions, type_collect)
                    type_ = suffix[0]
                    type_ = '境外输入' + type_ if abroad else type_
                    number = 0
                elif re.findall('在', phrase):
                    suffix = re.findall(r'(?<=在)[\D]+', phrase)
                    position = suffix[0]
                    if numbers:
                        # assert numbers[i - 1] == 1
                        number = numbers[i - 1]
                    else:
                        continue
                    type_ = types[i - 1]
                else:
                    continue
                positions.append(position)
                types.append(type_)
                numbers.append(number)
                i += 1
            result = pd.concat([result, pd.DataFrame([[type_collect] * len(types), types, positions, numbers],
                                                     index=['收集类型', '病例类型', '地点', '病例数/人数']).T])
    result.to_excel(f'{output_dir}/{date}.xlsx', index=False)
    print('Finish!')


if __name__ == '__main__':
    if sys.argv[1:]:
        try:
            extract_table(*sys.argv[1:])
        except Exception as e:
            print(f'{e}\nUsage:\n\tpython extract.py extract_text_file [save_dir]')
    else:
        read_text = 'test.txt'
        output_dir = './res'
        extract_table(read_text, output_dir)
