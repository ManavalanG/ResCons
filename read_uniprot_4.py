try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

hotspots_1 = ['metal ion-binding site', 'binding site', 'active site', 'modified residue', 'site', 'DNA-binding region']
hotspots_2 = ['sequence variant', 'mutagenesis site']
hotspots_all = hotspots_1 + hotspots_2
details = ['type', 'description']

out_handle = open('out.csv', 'w')
title_line = 'S.No,Primary Gene Name,Synonymous gene name,Starting aa,Ending aa,Hotspot type,Description,Original aa,Variant aa\n'
out_handle.write(title_line)

input_file = 'human/test.xml'
tree = ET.ElementTree(file=input_file)
root = tree.getroot()
entries = tree.findall('entry')

info_csv = ''
entry_no = 1
for entry in entries:
    primary_gene_name = None
    synonymous_gene_name = None

    for elem in entry.iterfind('gene/name'):
        if elem.attrib['type'] == 'primary':
            primary_gene_name = elem.text
            # print primary_gene_name
        elif elem.attrib['type'] == 'synonym':
            if synonymous_gene_name:
                synonymous_gene_name += ( ',' + elem.text)
            else:
                synonymous_gene_name = elem.text
        # else:
        #     primary_gene_name = elem.text
    print entry_no, primary_gene_name


    for elem in entry.iter():
        if elem.tag == 'feature':
            # print elem.tag, '1', elem.attrib, '2'
            checkpoint = 0
            temp = ''
            if elem.attrib['type'] in hotspots_all:
                for children in list(elem.iter()):
                    end_pos = False

                    for item in details:
                        if item in children.attrib:
                            temp += ',"%s"' %children.attrib[item]

                    if elem.attrib['type'] in hotspots_2:
                        if children.tag == 'original':
                            original = children.text
                        elif children.tag == 'variation':
                            variation = children.text


                    if children.tag == 'position':
                        begin_pos = int( children.attrib['position'] )
                        end_pos = int( children.attrib['position'] )
                    elif children.tag == 'begin':
                        begin_pos = int(children.attrib['position'])
                    elif children.tag == 'end':
                        end_pos = int(children.attrib['position'])

                    if end_pos:
                        if elem.attrib['type'] in hotspots_1:
                            info_csv = '%i,"%s","%s",%i,%i%s' %(entry_no, primary_gene_name, synonymous_gene_name, begin_pos, end_pos, temp)
                            out_handle.write(info_csv + '\n')

                        elif elem.attrib['type'] in hotspots_2:
                            info_csv = '%i,"%s","%s",%i,%i%s,"%s","%s"' %(entry_no, primary_gene_name, synonymous_gene_name, begin_pos, end_pos, temp, original, variation)
                            out_handle.write(info_csv + '\n')


    entry_no += 1


out_handle.close()