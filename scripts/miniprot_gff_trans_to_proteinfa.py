#!/usr/bin/env python

"""
Program: miniprot_gff_trans_to_proteinfa.py

This program takes a gff file from miniprot, output with the trans option, and exports
 the protein sequences in fasta format. The header of the protein sequence is the MP number,
 followed by a space, then the original protein sequence.

For this output:

##gff-version 3
##gff-version 3
##PAF   XP_033724493.1  368     119     368     +       OZ024763.1      30684678        14746944        14747688        498     783     0       AS:i:765        ms:i:765        np:i:193        fs:i:0  st:i:0  da:i:939        do:i:0  cg:Z:124M1I1M1I6M1I15M1D13M1I4M1D5M1I4M2I4M6I18M6D11M4D31M      cs:Z:*aacD:3*gatY*ctgS:2*aacD:4*tatF:13*agcG*gagD:12*gtcI*aacS:6*aagR:33*tacF:2*gttL:3*tacF*gaaN*ggtD:13*aaaR:4*acgA*acgA:3*cacP:6*tccN+H:1+G*accP*ttcY:4+Q:3*aacP:4*tccP*tacP*gcgY*cacS:3-cac:2*ccaS:2*acgA:1*aacT*ccgH*tccN*tccA*aacS:1+H*gccQ:3-cac:3*gtgI:1+N*tacF*agcN:2+TQ*cagA:3+PPGSTL:1*accS*ttcC:1*gccY*tacN:1*acaL*gccP*gccS*gccG*atgV:1*tcgL:1*atgP:1*ggcA-atgcgtatcgcccagccg*ccgG:1*aacT:1*ccaS:5*tcgN-ccggcgaccgcc*gccS*acgV:2*atgS*tcgG:2*gccP:1*accS:1*cagY*caaR:1*gccQ*gagG*ccgE:1*ttcI:2*atgV*cccH:4*gcaG:2
##STA   NSNKDLDPNQKPPYSYVALIAMAIKESSEKRLTLSGIYQYIVNKFPYYEKNKKGWQNSIRHNLSLNECFVKVPREGGGERKGNYWTVDPAYEGMFEKGNYRRRRRMKRPYRTTISLHKPLFADSCTFNQFALTKNYFSPSYAHQYSHQYPPWTLNPSSNAAGMGHMSQVSYSSCQRVPSTFGAYPTAAAMQSSMSGMRIAQPPTNYPQLNDYSPATAATSPMSPFAFTPQQQAEPSFNAMPYTYWAER*
OZ024763.1      miniprot        mRNA    14746945        14747691        765     +       .       ID=MP000001;Rank=1;Identity=0.6360;Positive=0.7395;Target=XP_033724493.1 120 368
OZ024763.1      miniprot        CDS     14746945        14747691        765     +       0       Parent=MP000001;Rank=1;Identity=0.6360;Target=XP_033724493.1 120 368
OZ024763.1      miniprot        stop_codon      14747689        14747691        0       +       0       Parent=MP000001;Rank=1
##PAF   XP_033724494.1  1325    19      1266    +       OZ024742.1      378904118       58210664        58214467        1356    3907    0       AS:i:1954       ms:i:1945       np:i:717        fs:i:4  st:i:0  da:i:60 do:i:111        cg:Z:50M2D8M1I100M3I13M1D9M5D13M4D9M3D34M2D13M17D47M2F34M1I36M2I35M1I314M4I57M1G5M8D5M1F47M3I46M2I38M2D2M1F51M1D114M4D25M1D18M10I5M7I35M2D20M1D26M1D3M cs:Z::1*aaaT*aagS:1*aaaF:1*gaaR:1*aatE*ctcR:1*ttgR:2*acaS:1*tgtL*gacL*ataT:1*ccaD*gacE*aaaV:1*aaaQ*actV*tcaN*atgT*tttL*ttaI*actY:1*atcM:1*tccG*aaaE*ggcA*cgaE*gacN:1*catV*tctK*actS*ctgF*gatG*ttcM*gatS*tcaA*cccA:1-gaagaa:1*aacK*cttK*caaY:2*atcV*gaaS+A*aaaR:1*gatN*gcaG:1*tgtF*gaaI*cccT:1*aaaR:1*ataT:2*ctgE:1*cacA:2*ttcN*acaH*tgtR*ggaK:1*gcgN:1*catG*caaE*aaaA*ttcA:1*gagS*tacF*gtgI*gtaT*gagA:1*cgcY*caaG*aagL*gcaV*caaE*cagH:2*ttcY:1*gacQ:1*acaR:1*agtE*ctcM*accI:2*attR*ctgI*ataV*tctV:4*atgS*agaK:1*agaS:1*agaK:1*ttgQ*agaL*gaaD*ccaA*aatT:1*gacT:1*caaE:1*acaA*attV*caaA*gcaQ*ggtV*cagR*aatQ*gcaS:1*caaA*acaV*cgaK*cgcK:1*tcaQ*aggT*atgV*ctgV*aatR+DSV*gccV:1*tccT*aaaT:1*agaP*gatS*aacS*gaaS*atcV*tcgD:1*ataL-cat:1*caaK*cacN*gacS*agaK:2*gccY:1-acacagataactcgc:1*atgH*caaG:1*agaK:2*aacY*gagD*caaK*ggtS*ggaH*agtT-gatagaaagccc*gcgS:3*gcgG*ccaQ*agtA*caaN*cacQ-tacgtcact*aacE:1*cacQ*tatR:2*tacF*aatK:1*caaR*agaK*aatD:3*acaA:1*cacN*aagA*ataT:1*gaaR:2*aagS:3*tacH*tatF*gcaE:1*atgK:1-aaatca:1*aaaT*tttT*gtaL*catS*aaaE*gtgI:1*acaQ*catY*tccK:2-gaatcagatgacgatatacatgcactcaaagtgatcactttaaataaatcc*cgcG*gtgY*catQ*aaaF*atcL*gaaG*aagS*ttcI*aagD*actD*acaS*aaaS*gacK*attP:3*gaaS:1*aaaC*gtgM*aatD*aacG*tcgL*aatP:1*acaE:10*gtaA*ctgI:1*ttcE*agaV*gaaT*tatF*aaaN*aatK*attL-at*tccK*aaaG*accL*ccaE*actL*gcaE:3*actK*cgaL:1*tcaQ*gccG*tacP*aatG*aacK*acaN*gctR*ataL*aaaR:4*tgcF*acaM*gcaS*gtcT*gttL*aaaT*tacT*aaaE:1+G:1*tcaQ*aaaA*caaT:1*aaaE*ttcV*attY:1*gcaI:2*aacL*tcaN*cacR:1*gttL*ataL*gatG*gcaR*aagP*acaA*gcaI*catE*gatK:3*attV*cagK*agaT*attV*ttaA*aaaD:1*aatH+QQ:1*tcgV*aacQ*gcaI*aacE*gacE*ataK*acaK*gaaK*aaaQ*tatF*gaaP*gacK*tgtL:1*ggaS*gaaG*atcL:1*actK*ctgF*ccaH*ataG*acaE*catY*catK:1*acaE*gtgL*gatR*ccaE*caaD*ataA*actK:1+Y*gtaA:1*catT*ccaT*ccaA:1*tcgR:1*cctA*tttL*tctP*ataL*aaaL*gaaD:1*ctaV:1*aaaD:2*aaaQ*aaaR:1*acaE*gaaQ*ttaQ:1*atcV:1*aaaS:2*agtV:5*gtaC*aatS*tccG:2*atcV:1*gaaP:1*gccK*aacD:1*aaaR*ctcV:1*ctcI:1*ttgV:1*cccL*aagT*aatQ:3*gccS*ataV*aaaR:1*caaE:1*ttaH*aaaQ:2*acaS*gcgV:1*gacQ*attT*tttL*gcaH*caaS*atgL*catT*aaaG*cctA:1*tatV:3*ctcI:2*tctK*aatA:1*tacF:4*gttL*gacE*gaaQ*gaaT:1*tccR*cacL:1*ttaT:2*aacI:6*cacC:1*acaN:3*tatF:2*cacS:2*agtP:1*gttH:2*aaaR*gaaR*atcM*gcaS*aaaQ*ataL*ataL:1*ggcS*ataE*gaaL:1*gcaV*agaV*aatC*tcaQ*cagM:2*attE*ataL:1*tggF*tccG*agcT*actN*ttaQ:3*acaD*aatQ:2*caaE*caaS:1*ttcL*gaaK*cgaK*atcL*cgaE*aagT*cacA:1*ttaV*aagT:2*aaaE*tcaE:2*atcE:1*aatS*gcaK*tcaP:1*ataV:1*tacF*cttV:2*ctaI*ttaI*accD:1*gaaT:2*aaaR*caaP:2*accK:5*attL*cacK*atgF*ccaE*ctaA:1*acaS*ggcN*aaaV*aaaS:1*ttaV*caaR:2*ctaM:1*atgI*acaV:1*tacQ:1*tgtG*aaaR:1*ttaS:2*tacI*tccA:1*gttK*acaS*gcaK:1*ttaI:1*cagD:3*agtK:2*ttaT:1*tcaL*ctcW*gacG*aaaP:2*atcS*gaaS:1*ataF*aaaQ*aagQ*ctaI:1*gaaQ*atgE*attL*accS*agaS*aaaT*ccaQ:2*caaS*tttL:3*gatA*ctaA:1*gtaT:1*ataV*actS*acaA:3*aagS*catY:3*gctG*ataV:1*gagL:2+PRDN*cacS*aacD*gaaQ*agaD:1*ctaK:2*tcgA*tttY*acaA:1*catR:1*acaM*actN*caaA*gcgT:2*aacR:1*tctA*ccgQ*ataV:1*agcK:3*agtA*atgA*gttT*ttcW:2*aagE:2*catS*cagD*tttL*ataL*tacV:1*actM:2*agaQ:2*aacT:3*gcgP:1*aV:1*aaaL*tatI*aaaG*aatS-aaaaaaataaaaaaaaaaaaaaaa:1*tcgP*atcL*actS*aaaE-t*gccL:3*atcV:3*ctaR*ctaM*gcaR:1*caaM:3*tttY:1*ctaI:1*tacH*ataV:3*gaaL*tcaL*gttY*gttI:2*actA:4*tacP:1*gaaR:1*aatS*actE:1*gaaP:1*tcaE*gacE+SFH*gaaQ*gatE*atgV*aaaE*tcgA*cacY*ataV*catD*gcaS*gtaI*ataL*aaaM:1*tacV*cacP*ataA*accS*gacE*actA*aaaR:1*aacE:1*tacI*aaaR*gaaL*gccK*actL*gcaR*cacE:2*tcaV*ctaC*cagG*acaI*ctcI*atcL*acaS:1*acaC*agaQ*acaD*gaaG:2+GY*ccaE:1*gacP*aaaS:1*ccaS*acaV*aacA*attT*aaaR:1*ttcY*gtcW*tctQ:1*cgtK:2*ataL:1*atcV*gtcC*aatQ:2*gttL:1*aaaR*agtG*ccaN:1*atcL*attV*gtaI:1*aagT:2-cataga*gaaR:1-a:1*gctI*ctaM:1*gaaK*ataL:1*accD:2*atgQ:2*gaaV:1*tacC*acaR:4*gaaD*tgtS*ataV*tacW:3*attM:1*gctR*caaE:1*aaaA:1*atgV*attV*aacK*acaN:1*tcaT*tacT:2*gatR*tacK*aaaR*aatA*caaD*caaI:1-gcc:3*ataK*aatP*catS:1*ataF:1*actE*acaR:2*accQ:1*gtaL:3*ataL:2*ctgW*tctK:1*aatS*ttgN:1*gtcL*acaL*atcV*attV:2*acaF*acaS:1*tttY*tttI*gatE:1*cacA*gaaK*ataL*aatS*gatS*tgcT*cagT:1*aaaP*actD*gtaI*atcV*aaaL*agaH*ataL*aaaQ*gaaS*actM:2*aaaR*ttcH:3*caaE:1*gttL:6*gaaQ*ttcY*aatA:2*gaaV:1*aaaQ*aaaQ:2*aaaD*aatK*tggY*catG:1*tatI:1*acgR*aaaT*acaS:3*tatH*catP:1*gcaG:2*atgE*gttA:2*aatA*ataV:2*ataV-aaaaagacgtta:1*aagR*gcaL:1*aaaG*agtG*ggcS*gaaK*gaaD*cctS*caaY:5*ttaY:1*acgA:1*aaaP:1*caaR*aatC:1-gca*ccaH:3*aatE*caaL*ataL:1*aatG:1*acaK*ttgI:2*aatT:2*aacE+NPNNLEPKWP*ataD:1*cacT:1*aaaF+RAKEREH*attR:2*aggQ*ccaK*cacQ*gaaN*cagF*aagD*aaaR:1*aaaQ*attR*acaA*acaK:1*ctcR:1*gagP:1*aaaR*ccaS*cacG:1*acgM:1*agaW*attV*aagS*ttcG*gatP*aacG*aatG*tggG*tttY-agaagg*ggcS*aaaG*ataT:1*aaaR*aaaE*ctaP*gaaV*caaS*cctD:2*tacF:2*gtaN*acaG*gagP*gaaS:1-aat:2*aaaR:3*caaR:1*attL*cttM*aaaL*acaE*aaaP:2*attL*cgcE*ataP*acaD:1*gaaP*ataE*gatG*tatM*gaaH:1-ata*cttF:2
##STA   WKKWKKEFNLFLIATECDIKPDKVKTSMFLTSIGSKGRDIHSTLDFDSPDEEMNLQTVIEKFDAYCEPRKNIIFLRHKFFTCGQAEHQKFDEYVVELRQKAQQCEFGDLTDSLTRDILISGIRDMRLRERLLREPNLDLQKTIQAGQNAEQTRRQSRMLNATSKTRDNEISAIHSQHDRKGAKTQITRHMQSRRQNEQGGSDRKPAKFQAPSQHYVTNCHYCGYNHQRNKCPTFHKICEKCKKKGYYAKMCKSRKFVHKVETHSDDESDDDIHALKVITLNKSRVHKIEKFKTTKDIWTVELKVNNSNVTFKIDTGADVTVLPFREYKNISKTPTAKTSTRLSAYNNTAIKVLGKCTAVVKYKKVSKQTKFIVAKTNSHAVIDAKTAHDLNLIQRILKINSSNANDITEKYEDCFGEIGTLPITHHITVDPQITPVIHPPRSVPFSIKEKLKKELKKMTELGIIKPVSEPTDWVNSLVIVEKANGKLRLCLDPKNLNKAIKRQHLKLPTAEDIFAQMHKPKYFSKLDASNGYWQIPVDEESSHLLTFNTPFGRYHFTRLPYGIHSASEVFQKEIAKIIDGIEGARNSQDDIIIWSSTLEEHTNRLQQVFERIRKHGLKLNKSKCIFNASQITYLGHLLTSEGIKQDPTKVQAIIHMPLPTGKKELQRFLGMTNYLCKFLPNYSEVTAPLRQLLQSDTLWSLDKPQIEAIKKLKEMITRKPVLQFYDPDLPVKITTDASKHGLGAILEQKHNERWLPVSFTSHATTQAEQNYSPIESEALSMVFACKKFHQFIYGTHFRIENDHKALSKYKNKKIKKKKKKSITKAPPRIQRFLLALQRYDFDLNYIPGKESVVADTLSRAYLEQNTPEVSDEDMKSHIHAVIKNYHITDTKLNEYKEATAHDSSLQTLITYTRTEWPPKDKIPTNIKPFVSVRDEISIVNGLVLKSPRIIVPKSLHREAEALQEIHTGHMGIEKYTERAKECIYWPGINAQIKDMINTCSYCIDYKNQQPAEPLINHEIPTTPWTKVGTDIFHLSGNLYVTIIDYTTRFFDLHEINDCQSKTVIKRIKETFAKFGIPQTVVSDNGPEFNSAEFKKFAKNWHFYHTKTSPHYHQANGMVERNIQTIKKTLKKAFKSGEEPQLALLALRTTKLQNGAPSPANQIMNRTLRTNLPNILHDKIVLRPHEQKKSKITTELPELKPHDTVRIKFDNNWFRRGKIIKKLEQPRSYLVVTEEGNTVKRNRQHILKTKENIRITEEIDYENILPNS
OZ024742.1      miniprot        mRNA    58210665        58214467        1945    +       .       ID=MP000002;Rank=1;Identity=0.3471;Positive=0.5506;Frameshift=4;Target=XP_033724494.1 20 1266
OZ024742.1      miniprot        CDS     58210665        58214467        1945    +       0       Parent=MP000002;Rank=1;Identity=0.3471;Frameshift=4;Target=XP_033724494.1 20 1266
##PAF   XP_033724497.1  276     160     194     -       OZ024740.1      442868002       338788896       338843969       39      102     0       AS:i:53 ms:i:88 np:i:23 fs:i:0  st:i:0  da:i:-1 do:i:60 cg:Z:19M54974U14M       cs:Z:*aagR:1*gttI*ggtN:1*aatQ:1*ggcD*aatH*agtA*gggA:1*ttcS:1*aatD*tatC:1*ggtA:1*gK~gt54971ag-ga*gttK*tctE*ccaR*cagE*attV:2*atcV:1*ttcY:3*ctcM
##STA   KCVGCNGGNSGSFRNYSGWGVSPQILRIKFTENLF
OZ024740.1      miniprot        mRNA    338788897       338843969       88      -       .       ID=MP000003;Rank=1;Identity=0.3824;Positive=0.6765;Target=XP_033724497.1 161 194
OZ024740.1      miniprot        CDS     338843912       338843969       49      -       0       Parent=MP000003;Rank=1;Identity=0.3500;Target=XP_033724497.1 161 179
OZ024740.1      miniprot        CDS     338788897       338788940       39      -       2       Parent=MP000003;Rank=1;Identity=0.4286;Target=XP_033724497.1 180 194

The resulting fasta entry will be:
>MP000001 XP_033724493.1
NSNKDLDPNQKPPYSYVALIAMAIKESSEKRLTLSGIYQYIVNKFPYYEKNKKGWQNSIR
HNLSLNECFVKVPREGGGERKGNYWTVDPAYEGMFEKGNYRRRRRMKRPYRTTISLHKPL
FADSCTFNQFALTKNYFSPSYAHQYSHQYPPWTLNPSSNAAGMGHMSQVSYSSCQRVPST
FGAYPTAAAMQSSMSGMRIAQPPTNYPQLNDYSPATAATSPMSPFAFTPQQQAEPSFNAM
PYTYWAER*
>MP000002 XP_033724494.1
WKKWKKEFNLFLIATECDIKPDKVKTSMFLTSIGSKGRDIHSTLDFDSPDEEMNLQTVIE
KFDAYCEPRKNIIFLRHKFFTCGQAEHQKFDEYVVELRQKAQQCEFGDLTDSLTRDILIS
GIRDMRLRERLLREPNLDLQKTIQAGQNAEQTRRQSRMLNATSKTRDNEISAIHSQHDRK
GAKTQITRHMQSRRQNEQGGSDRKPAKFQAPSQHYVTNCHYCGYNHQRNKCPTFHKICEK
CKKKGYYAKMCKSRKFVHKVETHSDDESDDDIHALKVITLNKSRVHKIEKFKTTKDIWTV
ELKVNNSNVTFKIDTGADVTVLPFREYKNISKTPTAKTSTRLSAYNNTAIKVLGKCTAVV
KYKKVSKQTKFIVAKTNSHAVIDAKTAHDLNLIQRILKINSSNANDITEKYEDCFGEIGT
LPITHHITVDPQITPVIHPPRSVPFSIKEKLKKELKKMTELGIIKPVSEPTDWVNSLVIV
EKANGKLRLCLDPKNLNKAIKRQHLKLPTAEDIFAQMHKPKYFSKLDASNGYWQIPVDEE
SSHLLTFNTPFGRYHFTRLPYGIHSASEVFQKEIAKIIDGIEGARNSQDDIIIWSSTLEE
HTNRLQQVFERIRKHGLKLNKSKCIFNASQITYLGHLLTSEGIKQDPTKVQAIIHMPLPT
GKKELQRFLGMTNYLCKFLPNYSEVTAPLRQLLQSDTLWSLDKPQIEAIKKLKEMITRKP
VLQFYDPDLPVKITTDASKHGLGAILEQKHNERWLPVSFTSHATTQAEQNYSPIESEALS
MVFACKKFHQFIYGTHFRIENDHKALSKYKNKKIKKKKKKSITKAPPRIQRFLLALQRYD
FDLNYIPGKESVVADTLSRAYLEQNTPEVSDEDMKSHIHAVIKNYHITDTKLNEYKEATA
HDSSLQTLITYTRTEWPPKDKIPTNIKPFVSVRDEISIVNGLVLKSPRIIVPKSLHREAE
ALQEIHTGHMGIEKYTERAKECIYWPGINAQIKDMINTCSYCIDYKNQQPAEPLINHEIP
TTPWTKVGTDIFHLSGNLYVTIIDYTTRFFDLHEINDCQSKTVIKRIKETFAKFGIPQTV
VSDNGPEFNSAEFKKFAKNWHFYHTKTSPHYHQANGMVERNIQTIKKTLKKAFKSGEEPQ
LALLALRTTKLQNGAPSPANQIMNRTLRTNLPNILHDKIVLRPHEQKKSKITTELPELKP
HDTVRIKFDNNWFRRGKIIKKLEQPRSYLVVTEEGNTVKRNRQHILKTKENIRITEEIDY
ENILPNS
>MP000003 XP_033724497.1
KCVGCNGGNSGSFRNYSGWGVSPQILRIKFTENLF
"""

import sys

def parse_gff(file_path):
    # Each entry of sequences is:
    #   {"proteinID": "protein_id_string",
    #    "sequence": "YASDKJHAUYGHKAJSHF*",
    #    "MP_id": "MP000001"}
    # Each entry forms one fasta entry
    sequences = []
    collecting_sequence = False

    current_entry = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('##PAF'):
                # We need to reset current_entry. If there is something in it, collect.
                if len(current_entry) == 3:
                    sequences.append(current_entry)
                # Use this to keep track of collecting.
                current_entry = {}
                parts = line.split('\t')
                if len(parts) > 1:
                    current_protein = parts[1]
                current_entry['proteinID'] = current_protein
                continue

            if line.startswith('##STA'):
                fields = line.split()
                if len(fields) != 2:
                    raise ValueError(f"Invalid sequence line: {line}")
                sequence = fields[1]
                current_entry["sequence"] = sequence
                continue

            # make sure that ID= is in the line
            if "ID=" in line:
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                attributes = {attr.split('=')[0]: attr.split('=')[1] for attr in fields[8].split(';') if '=' in attr}
                if "ID" in attributes:
                    current_entry["MP_id"] = attributes["ID"]

    # edge case to return the last one after the for loop
    if len(current_entry) == 3:
        sequences.append(current_entry)

    return sequences

def write_fasta(sequences, output_file):
    with open(output_file, 'w') as out:
        for entry in sequences:
            header = f">{entry['MP_id']} {entry['proteinID']}\n"
            out.write(header)
            for i in range(0, len(entry["sequence"]), 60):
                out.write(entry["sequence"][i:i+60] + "\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: miniprot_gff_trans_to_proteinfa.py <input.gff> <output.fasta>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    sequences = parse_gff(input_file)
    write_fasta(sequences, output_file)

if __name__ == '__main__':
    main()