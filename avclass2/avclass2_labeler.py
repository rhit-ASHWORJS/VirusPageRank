#!/usr/bin/env python3
'''
AVClass2 labeler
'''

import os
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1, os.path.join(script_dir, 'lib/'))
sys.path.insert(1, os.path.join(script_dir, '../shared/'))
import argparse
from avclass2_common import AvLabels
from operator import itemgetter
import evaluate_clustering as ec
import json
import traceback
import gzip

# Default tagging file
default_tag_file = os.path.join(script_dir, "data/default.tagging")
# Default expansion file
default_exp_file = os.path.join(script_dir, "data/default.expansion")
# Default taxonomy file
default_tax_file = os.path.join(script_dir, "data/default.taxonomy")

def guess_hash(h):
    ''' Given a hash string, guess the hash type based on the string length '''
    hlen = len(h)
    if hlen == 32:
        return 'md5'
    elif hlen == 40:
        return 'sha1'
    elif hlen == 64:
        return 'sha256'
    else:
        return None

def format_tag_pairs(l, taxonomy=None):
    ''' Return ranked tags as string '''
    if not l:
        return ""
    if taxonomy is not None:
        p = taxonomy.get_path(l[0][0])
    else:
        p = l[0][0]
    out = "%s|%d" % (p, l[0][1])
    for (t,s) in l[1:]:
        if taxonomy is not None:
            p = taxonomy.get_path(t) 
        else:
            p = t
        out += ",%s|%d" % (p, s)
    return out

def list_str(l, sep=", ", prefix=""):
    ''' Return list as a string '''
    if not l:
        return ""
    out = prefix + l[0]
    for s in l[1:]:
        out = out + sep + s
    return out

def main(args):
    # Select hash used to identify sample, by default MD5
    hash_type = args.hash if args.hash else 'md5'

    # If ground truth provided, read it from file
    gt_dict = {}
    if args.gt:
        with open(args.gt, 'r') as gt_fd:
            for line in gt_fd:
                gt_hash, family = map(str, line.strip().split('\t', 1))
                gt_dict[gt_hash] = family

        # Guess type of hash in ground truth file
        hash_type = guess_hash(list(gt_dict.keys())[0])

    # Create AvLabels object
    av_labels = AvLabels(args.tag, args.exp, args.tax,
                         args.av, args.aliasdetect)

    # Build list of input files
    # NOTE: duplicate input files are not removed
    ifile_l = []
    if (args.vt):
        ifile_l += args.vt
        ifile_are_vt = True
    if (args.lb):
        ifile_l += args.lb
        ifile_are_vt = False
    if (args.vtdir):
        ifile_l += [os.path.join(args.vtdir, 
                                  f) for f in os.listdir(args.vtdir)]
        ifile_are_vt = True
    if (args.lbdir):
        ifile_l += [os.path.join(args.lbdir, 
                                  f) for f in os.listdir(args.lbdir)]
        ifile_are_vt = False

    # Select correct sample info extraction function
    if not ifile_are_vt:
        get_sample_info = av_labels.get_sample_info_lb
    elif args.vt3:
        get_sample_info = av_labels.get_sample_info_vt_v3
    else:
        get_sample_info = av_labels.get_sample_info_vt_v2

    # Select output prefix
    out_prefix = os.path.basename(os.path.splitext(ifile_l[0])[0])

    # Initialize state
    first_token_dict = {}
    token_count_map = {}
    pair_count_map = {}
    vt_all = 0
    avtags_dict = {}
    stats = {'samples': 0, 'noscans': 0, 'tagged': 0, 'maltagged': 0,
             'FAM': 0, 'CLASS': 0, 'BEH': 0, 'FILE': 0, 'UNK': 0}

    # Process each input file
    for ifile in ifile_l:
        # Open file
        if args.gzip:
            fd = gzip.open(ifile, 'rt')
        else:
            fd = open(ifile, 'r')

        # Debug info, file processed
        sys.stderr.write('[-] Processing input file %s\n' % ifile)

        processed = 0
        viruses = []
        antiviruses = []

        votedict = {}

        # Process all lines in file
        for line in fd:

            # If blank line, skip
            if line == '\n':
                continue

            # Debug info
            if vt_all % 100 == 0:
                sys.stderr.write('\r[-] %d JSON read' % vt_all)
                sys.stderr.flush()
            vt_all += 1

            # Read JSON line
            vt_rep = json.loads(line)

            # Extract sample info
            sample_info = get_sample_info(vt_rep)

            # If no sample info, log error and continue
            if sample_info is None:
                try:
                    name = vt_rep['md5']
                    sys.stderr.write('\nNo scans for %s\n' % name)
                except KeyError:
                    sys.stderr.write('\nCould not process: %s\n' % line)
                sys.stderr.flush()
                stats['noscans'] += 1
                continue

            # Sample's name is selected hash type (md5 by default)
            name = getattr(sample_info, hash_type)

            # If the VT report has no AV labels, output and continue
            if not sample_info.labels:
                # sys.stdout.write('%s\t-\t[]\n' % (name))
                # sys.stderr.write('\nNo AV labels for %s\n' % name)
                # sys.stderr.flush()
                continue

            # Compute VT_Count
            vt_count = len(sample_info.labels)

            # Get the distinct tokens from all the av labels in the report
            # And print them. 
            try:
                av_tmp = av_labels.get_sample_tags(sample_info)
                tags = av_labels.rank_tags(av_tmp)
                

                # AV VENDORS PER TOKEN
                if args.avtags:
                    for t in av_tmp:
                        tmap = avtags_dict.get(t, {})
                        for av in av_tmp[t]:
                            ctr = tmap.get(av, 0)
                            tmap[av] = ctr + 1
                        avtags_dict[t] = tmap

                #Runnin
                if args.aliasdetect:
                    prev_tokens = set()
                    for entry in tags:
                        curr_tok = entry[0]
                        curr_count = token_count_map.get(curr_tok, 0)
                        token_count_map[curr_tok] = curr_count + 1
                        for prev_tok in prev_tokens:
                            if prev_tok < curr_tok:
                                pair = (prev_tok,curr_tok)
                            else:
                                pair = (curr_tok,prev_tok)
                            pair_count = pair_count_map.get(pair, 0)
                            pair_count_map[pair] = pair_count + 1
                        prev_tokens.add(curr_tok)

                # Collect stats
                # FIX: should iterate once over tags, 
                # for both stats and aliasdetect
                #runnin
                if tags:
                    stats["tagged"] += 1
                    if args.stats:
                        if (vt_count > 3):
                            stats["maltagged"] += 1
                            cat_map = {'FAM': False, 'CLASS': False,
                                       'BEH': False, 'FILE': False, 'UNK':
                                           False}
                            for t in tags:
                                path, cat = av_labels.taxonomy.get_info(t[0])
                                cat_map[cat] = True
                            for c in cat_map:
                                if cat_map[c]:
                                    stats[c] += 1

                # Check if sample is PUP, if requested
                if args.pup:
                    if av_labels.is_pup(tags, av_labels.taxonomy):
                        is_pup_str = "\t1"
                    else:
                        is_pup_str = "\t0"
                else:
                    is_pup_str =  ""

                # Select family for sample if needed,
                # i.e., for compatibility mode or for ground truth
                if args.c or args.gt:
                    fam = "SINGLETON:" + name
                    # fam = ''
                    for (t,s) in tags:
                        cat = av_labels.taxonomy.get_category(t)
                        if (cat == "UNK") or (cat == "FAM"):
                            fam = t
                            break

                # Get ground truth family, if available
                if args.gt:
                    first_token_dict[name] = fam
                    gt_family = '\t' + gt_dict.get(name, "")
                else:
                    gt_family = ""

                # Get VT tags as string
                if args.vtt:
                    vtt = list_str(sample_info.vt_tags, prefix="\t")
                else:
                    vtt = ""

                # Print family (and ground truth if available) to stdout
                #myprint
                if not args.c:
                    if args.path:
                        tag_str = format_tag_pairs(tags, av_labels.taxonomy)
                    else:
                        tag_str = format_tag_pairs(tags)
                    # sys.stdout.write('%s\t%d\t%s%s%s%s\n' %
                    #                  (name, vt_count, tag_str, gt_family,
                    #                   is_pup_str, vtt))
                    # print(tag_str)

                    if name not in viruses:
                        viruses.append(name)

                    for key in av_tmp.keys():
                        avlist = av_tmp[key]
                        for av in avlist:
                            if av not in antiviruses:
                                antiviruses.append(av)
                    
                    for tag in av_tmp.keys():
                        if tag in ['01c02300fa', '15872a', '180solutions', '1clickdownload', '1gwaasep03e', '2345explorer', '383media', '4fd2e94b', '4shared', '5616bd', '84f6ca', 'aabe', 'aahe', 'abcu', 'accphish', 'acvt', 'acwbdswuqlp', 'adaebook', 'adagent', 'adclicer', 'addlyrics', 'addrop', 'adduser', 'adgazelle', 'adhelper', 'adinstaller', 'adkor', 'adloader', 'adond', 'adsuproot', 'adsvc', 'adultbrowser', 'advancedpccare', 'advheur', 'adwapper', 'adwareprofessional', 'adwin', 'afio', 'agentb', 'agobot', 'airinstaller', 'aksula', 'aliser', 'allaple', 'alman', 'alnaddy', 'alphabet', 'aluigi', 'alyak', 'amonetize', 'antavmu', 'antilam', 'antivirus', 'antivirusapp', 'antivirushijack', 'anyprotect', 'appshake', 'archsms', 'ardamax', 'ardurk', 'aspxor', 'asterisk', 'atraps', 'atros', 'au0ba0o9aili', 'audhi', 'autinject', 'autit', 'avdown', 'avinfoassist', 'b7c2a', 'babonock', 'babylon', 'badjoke', 'badur', 'bagle', 'baidu', 'baidusearch', 'balisdat', 'bamgadin', 'bamital', 'banayu', 'banbra', 'bancos', 'bandoo', 'banload', 'banz', 'barys', 'baryumka', 'basun', 'bawswerps', 'bayrob', 'bbbn', 'bdgk', 'beastdoor', 'beaugrit', 'bebo', 'beebone', 'begman', 'begseabug', 'behav', 'benjamin', 'berbew', 'bertle', 'bestafera', 'betabot', 'betload', 'betterinternet', 'bettersurf', 'bezigate', 'bhouninstaller', 'bibei', 'bicololo', 'bifrose', 'bifrost', 'binder', 'bindex', 'bionet', 'bitcro', 'bitman', 'black', 'blackhole', 'blackshades', 'bladabindi', 'blafixel', 'blakamba', 'blanajog', 'blazgel', 'blihan', 'blohi', 'bloored', 'blueh', 'bluesoft', 'bmkexmzzylgi', 'bmmedia', 'boaxxe', 'bobax', 'bobrowser', 'boht', 'bomgar', 'boopcel', 'boran', 'bototer', 'bprotector', 'braban', 'brappware', 'bredavi', 'bredolab', 'brontok', 'browsec', 'browsecx', 'browsefox', 'browserguard', 'browserpassview', 'browshot', 'brrowho', 'bruteforce', 'btmagnat', 'bublik', 'bulta', 'bundleloader', 'bundlore', 'bunitu', 'buterat', 'buzus', 'buzy', 'bzub', 'c4dlmedia', 'cain', 'calac', 'calper', 'capredeam', 'capsfin', 'carberp', 'cardspy', 'cashfiesta', 'cassiopeia', 'catalina', 'ceeinject', 'cekar', 'celofot', 'centrumloader', 'cerber', 'cgau', 'cgub', 'cheatengine', 'chifrax', 'chinad', 'chinbo', 'chindo', 'chir', 'chisburg', 'chiton', 'cidox', 'cinmus', 'cleaman', 'clikug', 'clons', 'cloudguard', 'cmqabctwk6ci', 'cnbtech', 'cngp', 'cnnic', 'cnzzbot', 'codewall', 'codiby', 'cognosads', 'coidung', 'collected', 'colooader', 'comame', 'compete', 'computrace', 'comrerop', 'conduit', 'confuser', 'consmiper', 'constructor', 'convertad', 'coolmirage', 'cosmicduke', 'cosmu', 'cosne', 'cospet', 'cossta', 'coulomb', 'couponmarvel', 'coupons', 'cpete', 'cpush', 'craagle', 'crawler', 'crime', 'critloki', 'crossrider', 'crowti', 'cryakl', 'cryfile', 'crypmodadv', 'cryptolocker', 'cryptos', 'crytex', 'csdimonetize', 'cutwail', 'cwqh', 'cybergate', 'cycbot', 'dadobra', 'dagava', 'dalexis', 'danginex', 'dapato', 'darkkomet', 'dartsmound', 'dater', 'davobevix', 'daws', 'dealply', 'deepsea', 'defaulttab', 'delbar', 'delf', 'delfi', 'delfiles', 'delfinject', 'delflash', 'delphi', 'delud', 'demp', 'dervec', 'deshacop', 'detnat', 'detroie', 'dewnad', 'dfsm', 'dial', 'dialupass', 'diamin', 'digitala', 'dinwod', 'dipad', 'diplugem', 'directdownloader', 'disfa', 'diztakun', 'dlboost', 'dlhelper', 'dllinject', 'dluca', 'dnguard', 'dnschanger', 'dnsunlocker', 'dodiw', 'doedlid', 'dofoil', 'dollarrevenue', 'domaiq', 'dorando', 'dorf', 'dorgam', 'dorifel', 'dorkbot', 'dorn', 'dorv', 'dotdo', 'dowgav', 'downloadadmin', 'downloadassistant', 'downloadguide', 'downloadhelper', 'downloadsponsor', 'downvision', 'dprotect', 'dridex', 'driverupdate', 'drmsoft', 'drolnux', 'droma', 'dsbot', 'dudra', 'dumaru', 'dumpy', 'dunsenr', 'dupator', 'dvhf', 'dwis', 'dybalom', 'dycler', 'dynamer', 'dyre', 'easyspeedcheck', 'ebgx', 'edeals', 'edownloader', 'eggnog', 'egroupdial', 'eibdbsjvwn', 'eicartest', 'elex', 'elitemax', 'elzob', 'emager', 'emotet', 'enchanim', 'enigma', 'enigmaprotector', 'enistery', 'enterok', 'eorezo', 'epack', 'epicscale', 'ertfor', 'esfury', 'eszjuxuan', 'etranslatorpro', 'eupuds', 'excrevie', 'exent', 'exetemp', 'expiro', 'express', 'expressdownloader', 'expressinstaller', 'extbro', 'extcrome', 'extenbro', 'eyestye', 'eziriznetreactor', 'ezula', 'fakeal', 'fakefolder', 'fakeie', 'fakeinit', 'fakeinstaller', 'fakeqq', 'fakerean', 'fakespypro', 'fakesysdef', 'faketool', 'fakewindow', 'fareit', 'farfli', 'fasong', 'fayu', 'fearso', 'fednu', 'fenomengame', 'fesber', 'filefinder', 'filemonster', 'filetour', 'filoskeed', 'firseria', 'firstfloor', 'firstinj', 'fleercivet', 'flowspirit', 'flystudio', 'fmjfas7dkyc', 'folcom', 'forcestartpage', 'fosniw', 'fostsy', 'fourthrem', 'fraudload', 'fraudpack', 'fraudrop', 'fraudster', 'frelex', 'frethog', 'fsysna', 'fuck', 'fujacks', 'fulfillingapps', 'fullscreen', 'funapps', 'funlove', 'fusing', 'fusioncore', 'gamarue', 'game', 'gamehack', 'gamehuck', 'gamemodding', 'gametea', 'gamethief', 'gametool', 'gamevance', 'gamup', 'ganelp', 'ganipin', 'garrun', 'gatak', 'gator', 'gbdialer', 'gemius', 'gena', 'gendal', 'genkryptik', 'genocide', 'gepys', 'geral', 'gesong', 'gessto', 'getnow', 'gibi', 'gigaclicks', 'gimemo', 'gippers', 'girlinred', 'globalupdate', 'glox', 'glupteba', 'gmunpackerinstaller', 'gobot', 'goforfiles', 'gofot', 'goldun', 'golroted', 'goobzo', 'googupdate', 'goredir', 'goriadu', 'gorillaprice', 'gotango', 'grandmedia', 'graybird', 'grenam', 'guagua', 'gudra', 'hacdef', 'hackav', 'hafen', 'hakaglan', 'hamweq', 'harnig', 'haxspy', 'hbpw', 'hebogo', 'heim', 'heri', 'heye', 'hezhi', 'hfyi', 'hicrazyk', 'hiddencamera', 'hiddeninstall', 'hidebaid', 'hideexec', 'hideicon', 'hideproc', 'hidewindows', 'hijacker', 'hiloti', 'hlubea', 'hlux', 'hmgfasxx', 'hmtoolbar', 'hoax', 'homa', 'hometab', 'hongdawanfang', 'hopadef', 'horst', 'hosts', 'hotbar', 'hotdownloads', 'hpdefender', 'hrup', 'huhk', 'hupigon', 'hype', 'i1sk', 'ibryte', 'iciko', 'icloader', 'idlekms', 'ilcrypt', 'ilheur', 'iluspy', 'imali', 'imestartup', 'iminent', 'inbox', 'induc', 'innomod', 'inservice', 'installall', 'installbrain', 'installcore', 'installerex', 'installflash', 'installiq', 'installmate', 'installmetrix', 'installmonetizer', 'installmonster', 'installtoolbar', 'instally', 'invader', 'iobit', 'ipamor', 'ircbot', 'iroffer', 'irplan', 'istartsurf', 'iterign', 'itorrent', 'iwin', 'jaik', 'jaku', 'jakuz', 'jatif', 'jawego', 'jeefo', 'jevafus', 'jewdo', 'jintor', 'joiner', 'jongiti', 'joosoft', 'jottix', 'jrat', 'judo', 'juntador', 'jyfi', 'kapser', 'kasidet', 'kasinst', 'katar', 'kedebe', 'kelihos', 'keybase', 'keyloggeronline', 'keymake', 'keystart', 'kflate', 'kiayksayren', 'kiction', 'kilim', 'killproc', 'klez', 'klovbot', 'koblu', 'kolab', 'kolabc', 'komodia', 'koobface', 'koutodoor', 'kovter', 'kraddare', 'kranet', 'krunchy', 'kuaiba', 'kuluoz', 'kuping', 'laban', 'lamer', 'laqma', 'ldpinch', 'lebag', 'lenovo', 'lethic', 'liha', 'liimpact', 'limitail', 'lineage', 'linkoptimizer', 'linkular', 'linkury', 'lipler', 'llac', 'lmir', 'lnkwinkap', 'loader', 'loadmoney', 'loadwar', 'locky', 'lohocla', 'lola', 'lollipop', 'loring', 'loudmo', 'loveletter', 'lovesong', 'lowzones', 'ltlogger', 'luder', 'luhe', 'lunam', 'lydra', 'lynx', 'lyposit', 'lyrics', 'mabezat', 'madon', 'maener', 'magania', 'magicbit', 'mailpassview', 'malex', 'malwareromovalbot', 'malwares', 'mamianune', 'manbat', 'mantal', 'maozhi', 'mapdimp', 'mapler', 'mask', 'masta', 'matsnu', 'maxdriver', 'maxiget', 'maximus', 'mayachok', 'medfos', 'mediadrug', 'mediafinder', 'mediaget', 'mediamagnet', 'megasearch', 'meinhudong', 'memery', 'menti', 'mepaow', 'messengerplus', 'messo', 'mewsspy', 'mh0aawheovm', 'michela', 'microfake', 'microjoin', 'midia', 'midie', 'mikey', 'mimail', 'mimikatz', 'minari', 'minggy', 'miniduke', 'minilash', 'mira', 'misskaz', 'mizenota', 'mkar', 'mobler', 'mofin', 'mogoogwi', 'mokes', 'moleboxultra', 'moleboxvs', 'monder', 'monitoringtool', 'moonlight', 'mozla', 'mpress', 'msildrop', 'msilinj', 'msilkrypt', 'msilperseus', 'msposer', 'mudrop', 'multibar', 'multidl', 'multidropper', 'multipacked', 'multiplug', 'multitoolbar', 'murofet', 'musin', 'mutabaha', 'mydoom', 'mypcbackup', 'mytob', 'mytonel', 'mywebsearch', 'myxah', 'nagoot', 'nagram', 'nakuru', 'nanobot', 'nanocore', 'narcha', 'navipromo', 'nayrabot', 'nbdd', 'necast', 'necurs', 'neobar', 'neojit', 'nepoe', 'neshta', 'netcat', 'netins', 'netpass', 'netseal', 'netshrink', 'netsky', 'nettool', 'netwiredrc', 'neurevt', 'newloader', 'newplayer', 'nextlive', 'ngrbot', 'nieguide', 'nitol', 'nivasteg', 'nivdort', 'njrat', 'noadware', 'noancooe', 'noobyprotect', 'noplemento', 'nosok', 'nsanti', 'nsismod', 'nspm', 'ntrootkit', 'nuprader', 'nurech', 'nurjax', 'nymaim', 'observer', 'oc1bzuzfpooj', 'ociyota', 'ocna', 'ocsbundle', 'offend', 'offtoup', 'oficla', 'oleloa', 'olufus', 'omaneat', 'oneclickdownloader', 'onescan', 'onesystemcare', 'onlgame', 'onlinegames', 'opanki', 'opencandy', 'openinstall', 'opiker', 'optimizer', 'optimizerelitemax', 'optix', 'orbus', 'orcusrot', 'orsam', 'orto', 'otezinu', 'outbrowse', 'oxypumper', 'packedent', 'pakes', 'palevo', 'papras', 'parite', 'paskod', 'passfox', 'pasta', 'pastaleads', 'pcast', 'pcclient', 'pcfixcleaner', 'pchealthboost', 'pckeeper', 'pcmega', 'peermarket', 'pef13c', 'pennybee', 'penzievs', 'perflogger', 'perinet', 'perion', 'perkesh', 'personalsheild', 'pfoenic', 'pgpme', 'pher', 'phires', 'phorpiex', 'piccolor', 'picsys', 'pincav', 'pinguide', 'pintu', 'pioneer', 'piptea', 'pirminay', 'pistolar', 'pitit', 'pjtbinder', 'playmp3z', 'plemood', 'plimrost', 'plorexie', 'plsex', 'pluginaccess', 'plugx', 'podnuha', 'poison', 'polip', 'polycrypt', 'pondfull', 'popad', 'popdeals', 'pophot', 'pornoblocker', 'porntool', 'powp', 'pricegong', 'prifou', 'primecasino', 'privateexeprotector', 'processpatcher', 'procpatcher', 'proinstall', 'prorat', 'prosti', 'protector', 'protux', 'proxychanger', 'pskill', 'psoriasis', 
                        'pswtool', 'pugeju', 'pullupdate', 'puma', 'purityscan', 'pvlognetprotector', 'pwdump', 'pwsime', 'pwszbot', 'pykspa', 'qaccel', 'qadars', 'qbot', 'qeds', 'qhost', 'qihoo', 'qiwmonk', 'qjwmonkey', 'qqfish', 'qqgetpass', 'qqhelper', 'qqpass', 'qqrob', 'qqware', 'qsii', 'quervar', 'quireap', 'rack', 'radmin', 'radonskra', 'ramdo', 'ramnit', 'ranapama', 'ranbyus', 'ranos', 'ranserkd', 'ravs', 'rayra', 'razy', 'rbot', 'rebhip', 'recal', 'recam', 'recodrop', 'reconyc', 'recslurp', 'redaptor', 'redgirl', 'redosdru', 'redyms', 'refroso', 'registrybooster', 'regprocleaner', 'regrun', 'regsup', 'reimagerepair', 'reklosoft', 'relevantknowledge', 'relnek', 'remoteadmin', 'renaz', 'renos', 'resdro', 'resur', 'rewriteboota', 'ridnu', 'rimod', 'ripinip', 'rjump', 'rkproc', 'rlsloup', 'rockettab', 'rodecap', 'rofin', 'rombrast', 'rontokbro', 'ropest', 'rottentu', 'rovnix', 'rozena', 'rubar', 'rubinurd', 'rukap', 'rukoma', 'rukometa', 'runner', 'ruskill', 'rustock', 'sacto', 'sadenav', 'sahat', 'sality', 'samon', 'sanctionedmedia', 'sasfis', 'sasquor', 'sbwatchman', 'sbyinying', 'scar', 'scarpnex', 'scarsi', 'sciagnij', 'scramblewrapper', 'sdbot', 'sdld', 'searchaid', 'securityxploded', 'sefnit', 'sekur', 'selfish', 'senta', 'seodec', 'sfdld', 'sfone', 'sfuzuan', 'sfvwg', 'shade', 'shakblades', 'shaosmine', 'sharik', 'shark', 'shelma', 'shenzhen', 'shipup', 'shiz', 'shodi', 'shopper', 'shopperpro', 'shopperz', 'shopro', 'shouqu', 'shutdowner', 'shyape', 'shylock', 'siggen3ent', 'sikeaqnqmgjb', 'silcon', 'sileco', 'silly', 'sillyfdc', 'sillyw', 'simbot', 'simda', 'similagro', 'simplefiles', 'simplytech', 'sinowal', 'sinresby', 'sisbot', 'sisproc', 'sisron', 'sivis', 'sixer', 'skeeyah', 'skillis', 'skintrim', 'skor', 'sksocket', 'skyli', 'slenfbot', 'slimware', 'slmvwk', 'slugin', 'slym', 'smalltro', 'smartapps', 'smartassembly', 'smartfixer', 'smartfortress', 'smartinstaller', 'smee', 'smileonline', 'smshoax', 'smwnd', 'sober', 'socks', 'soft', 'soft2cn', 'soft32downloader', 'softcnapp', 'softobase', 'softomate', 'softonic', 'softpulse', 'sogou', 'sohanad', 'sohand', 'sohuva', 'somoto', 'sopinar', 'soulclose', 'sourtoff', 'soxgrave', 'spafx', 'spammy', 'speedbit', 'speedingupmypc', 'speil', 'spigot', 'spnr', 'sprotector', 'spybot', 'spyeye', 'spynet', 'spywarebot', 'squarenet', 'sram', 'staget', 'stagol', 'startsurf', 'staser', 'steam', 'stegvob', 'stimilik', 'stresid', 'stressreducer', 'strumapine', 'stub', 'subtab', 'superfish', 'supertuneup', 'surveyer', 'susnn', 'suweezy', 'svcminer', 'swisyn', 'swizzor', 'swrort', 'syddld', 'syncopate', 'sysn', 'systemcall', 'systemdefender', 'systemhealer', 'systemtool', 'systemtweaker', 'systex', 'systroj', 'systweak', 'sytro', 'tank', 'taobao', 'taojinstar', 'taranis', 'tc3ahj7nvlei', 'tcpscan', 'tdss', 'techsnab', 'temonde', 'tempedreve', 'temr', 'tencent', 'tenga', 'terkcop', 'teront', 'teslacrypt', 'testing', 'tibia', 'tibs', 'tinba', 'tirrip', 'tiyy', 'tmhfbm4amkj', 'tobfy', 'tofsee', 'toga', 'toggle', 'tonmye', 'topa', 'topmedia', 'toptools', 'torntv', 'torr', 'torrentclient', 'torrentsearch', 'toxic', 'tpyn', 'trojannotifier', 'truedownloader', 'trustezeb', 'trymedia', 'tufik', 'turbobit', 'turist', 'turkojan', 'tuscas', 'twexag', 'u07aby', 'ubibila', 'ubot', 'ugmk', 'ulpm', 'ultimatedefender', 'umnfrvwrr2pib', 'uniblue', 'unruy', 'upatre', 'urausy', 'urelas', 'ursnif', 'usteal', 'valcaryx', 'vanbot', 'vapsup', 'vasor', 'vbates', 'vbclone', 'vbinder', 'vbklog', 'vbkryjetor', 'vbran', 'vburses', 'vdld', 'vehidis', 'venik', 'verind', 'vermid', 'vernet', 'verti', 'veryhighconfidence', 'vidro', 'vietkey', 'viking', 'vilsel', 'vindor', 'virledi', 'virlock', 'virtumonde', 'virut', 'visicom', 'vitruvian', 'vittalia', 'vkhost', 'vkont', 'vmpbad', 'vmprotbad', 'vobfus', 'vonteera', 'vopackage', 'vopak', 'vorloma', 'vprotect', 'vtflooder', 'vundo', 'wabot', 'wajam', 'waldek', 'walta', 'wammuras', 'wapomi', 'warmup', 'watchman', 'wavipeg', 'wdjange', 'wdjiange', 'webalta', 'webbar', 'webcake', 'webdevaz', 'webfilter', 'webhat', 'webprefix', 'webtoos', 'webwatcher', 'wecod', 'weiduan', 'weird', 'wenper', 'wews', 'whiteice', 'widdit', 'widgi', 'wigon', 'wimg', 'winactivator', 'windef', 'winemm', 'winfetcher', 'wininf', 'winloadsda', 'winlock', 'winner', 'winsecsrv', 'wintrim', 'winvnc', 'winwebsec', 'winwrapper', 'wizzcaster', 'wizzlab', 'wlksm', 'wlord', 'wonton', 'wootbot', 'wowlik', 'wpakill', 'wprotmanager', 'writos', 'wsgame', 'wslhq', 'wswhacker', 'wuji', 'xadupi', 'xanfpezes', 'xetapp', 'xhnl', 'xiaohao', 'xiaoho', 'xiaoxiong', 'xiazai', 'ximera', 'xiquitir', 'xmkfrtraoceb', 'xorala', 'xorer', 'xorist', 'xpaj', 'xpyn', 'xtrat', 'xyligan', 'yabector', 'yabinder', 'yahlover', 'yakes', 'yantai', 'yappyz', 'yazz', 'yelloader', 'yellsob', 'yemrok', 'ymeta', 'yobdam', 'yoddos', 'yogosojo', 'yourinstaller', 'youxun', 'yqwaaasdfofb', 'yu0bfdsueddi', 'yuner', 'zafi', 'zanoza', 'zapchast', 'zbot', 'zboter', 'zeeborot', 'zegost', 'zemot', 'zepfod', 'zerber', 'zeroaccess', 'zestyfind', 'zhangguojian', 'zlob', 'zmlfa4s3unlb', 'zombienews', 'zonebac', 'zorg', 'zortob', 'zugo', 'zurgop', 'zusy', 'zuten', 'zvuzona', 'zxrloader', 'zygug', 'zylom', 'zzinfor']:
                            avlist = av_tmp[key]
                            for antivirus in avlist:
                                votedict[str(antivirus) + ":" + str(name)] = tag
                            # votedict["abc:123"] = "cool"
                    # print("processed: " + str(processed))
                    processed = processed + 1
                    # print(av_labels.taxonomy)

                # else:
                    # sys.stdout.write('%s\t%s%s%s\n' %
                    #                  (name, fam, gt_family, is_pup_str))
            except:
                traceback.print_exc(file=sys.stderr)
                continue
        
        sys.stdout.write("TopLeft")
        for virus in viruses:
            sys.stdout.write(","+virus)
        sys.stdout.write("\n")

        for antivirus in antiviruses:
            # for char in antivirus:
            #     sys.stdout.write(char)
            sys.stdout.write(antivirus.replace("—","-").replace("–","-"))
            for virus in viruses:
                if (antivirus + ":" + virus) in votedict.keys():
                    sys.stdout.write("," + votedict[antivirus + ":" + virus])
                else:
                    sys.stdout.write("," + "_")
            sys.stdout.write("\n")
        # Debug info
        # sys.stderr.write('\r[-] %d JSON read' % vt_all)
        # sys.stderr.flush()
        # sys.stderr.write('\n')

        # Close file
        fd.close()

    # Print statistics
    # sys.stderr.write(
    #         "[-] Samples: %d NoScans: %d NoTags: %d GroundTruth: %d\n" % (
    #             vt_all, stats['noscans'], vt_all - stats['tagged'], 
    #             len(gt_dict)))

    # If ground truth, print precision, recall, and F1-measure
    if args.gt:
        precision, recall, fmeasure = \
                    ec.eval_precision_recall_fmeasure(gt_dict,
                                                      first_token_dict)
        sys.stderr.write(
            "Precision: %.2f\tRecall: %.2f\tF1-Measure: %.2f\n" % \
                          (precision, recall, fmeasure))

    # Output stats
    if args.stats:
        stats_fd = open("%s.stats" % out_prefix, 'w')
        num_samples = vt_all
        stats_fd.write('Samples: %d\n' % num_samples)
        num_tagged = stats['tagged']
        frac = float(num_tagged) / float(num_samples) * 100
        stats_fd.write('Tagged (all): %d (%.01f%%)\n' % (num_tagged, frac))
        num_maltagged = stats['maltagged']
        frac = float(num_maltagged) / float(num_samples) * 100
        stats_fd.write('Tagged (VT>3): %d (%.01f%%)\n' % (num_maltagged, frac))
        for c in ['FILE','CLASS','BEH','FAM','UNK']:
            count = stats[c]
            frac = float(count) / float(num_maltagged) * 100
            stats_fd.write('%s: %d (%.01f%%)\n' % (c, stats[c], frac))
        stats_fd.close()

    # Output vendor info
    if args.avtags:
        avtags_fd = open("%s.avtags" % out_prefix, 'w')
        for t in sorted(avtags_dict.keys()):
            avtags_fd.write('%s\t' % t)
            pairs = sorted(avtags_dict[t].items(),
                            key=lambda pair : pair[1],
                            reverse=True)
            for pair in pairs:
                avtags_fd.write('%s|%d,' % (pair[0], pair[1]))
            avtags_fd.write('\n')
        avtags_fd.close()

    # If alias detection, print map
    if args.aliasdetect:
        # Open alias file
        alias_filename = out_prefix + '.alias'
        alias_fd = open(alias_filename, 'w+')
        # Sort token pairs by number of times they appear together
        sorted_pairs = sorted(
            pair_count_map.items(), key=itemgetter(1))
        # sorted_pairs = sorted(
        #     pair_count_map.items())

        # Output header line
        alias_fd.write("# t1\tt2\t|t1|\t|t2|\t"
                       "|t1^t2|\t|t1^t2|/|t1|\t|t1^t2|/|t2|\n")
        # Compute token pair statistic and output to alias file
        for (t1, t2), c in sorted_pairs:
            n1 = token_count_map[t1]
            n2 = token_count_map[t2]
            if (n1 < n2):
                x = t1
                y = t2
                xn = n1
                yn = n2
            else:
                x = t2
                y = t1
                xn = n2
                yn = n1
            f = float(c) / float(xn)
            finv = float(c) / float(yn)
            if args.path:
                x = av_labels.taxonomy.get_path(x)
                y = av_labels.taxonomy.get_path(y)
            alias_fd.write("%s\t%s\t%d\t%d\t%d\t%0.2f\t%0.2f\n" % (
                x, y, xn, yn, c, f, finv))
        # Close alias file
        alias_fd.close()
        sys.stderr.write('[-] Alias data in %s\n' % (alias_filename))


if __name__=='__main__':
    argparser = argparse.ArgumentParser(prog='avclass2_labeler',
        description='''Extracts tags for a set of samples.
            Also calculates precision and recall if ground truth available''')

    argparser.add_argument('-vt', action='append',
        help='file with VT reports '
             '(Can be provided multiple times)')

    argparser.add_argument('-lb', action='append',
        help='file with simplified JSON reports'
             '{md5,sha1,sha256,scan_date,av_labels} '
             '(Can be provided multiple times)')

    argparser.add_argument('-vtdir',
        help='existing directory with VT reports')

    argparser.add_argument('-lbdir',
        help='existing directory with simplified JSON reports')

    argparser.add_argument('-vt3', action='store_true',
        help='input are VT v3 files')

    argparser.add_argument('-gz', '--gzip',
        help='file with JSON reports is gzipped',
        action='store_true')

    argparser.add_argument('-gt',
        help='file with ground truth. '
             'If provided it evaluates clustering accuracy. '
             'Prints precision, recall, F1-measure.')

    argparser.add_argument('-vtt',
        help='Include VT tags in the output.',
        action='store_true')

    argparser.add_argument('-tag',
        help='file with tagging rules.',
        default = default_tag_file)

    argparser.add_argument('-tax',
        help='file with taxonomy.',
        default = default_tax_file)

    argparser.add_argument('-exp',
        help='file with expansion rules.',
        default = default_exp_file)

    argparser.add_argument('-av',
        help='file with list of AVs to use')

    argparser.add_argument('-avtags',
        help='extracts tags per av vendor',
        action='store_true')

    argparser.add_argument('-pup',
        action='store_true',
        help='if used each sample is classified as PUP or not')

    argparser.add_argument('-p', '--path',
        help='output.full path for tags',
        action='store_true')

    argparser.add_argument('-hash',
        help='hash used to name samples. Should match ground truth',
        choices=['md5', 'sha1', 'sha256'])

    argparser.add_argument('-c',
        help='Compatibility mode. Outputs results in AVClass format.',
        action='store_true')

    argparser.add_argument('-aliasdetect',
        action='store_true',
        help='if used produce aliases file at end')

    argparser.add_argument('-stats',
                           action='store_true',
                           help='if used produce 1 file '
                                'with stats per category '
                                '(File, Class, '
                                'Behavior, Family, Unclassified)')

    args = argparser.parse_args()

    if not args.vt and not args.lb and not args.vtdir and not args.lbdir:
        sys.stderr.write('One of the following 4 arguments is required: '
                          '-vt,-lb,-vtdir,-lbdir\n')
        exit(1)

    if (args.vt or args.vtdir) and (args.lb or args.lbdir):
        sys.stderr.write('Use either -vt/-vtdir or -lb/-lbdir. '
                          'Both types of input files cannot be combined.\n')
        exit(1)

    if args.tag:
        if args.tag == '/dev/null':
            sys.stderr.write('[-] Using no tagging rules\n')
        else:
            sys.stderr.write('[-] Using tagging rules in %s\n' % (
                              args.tag))
    else:
        sys.stderr.write('[-] Using default tagging rules in %s\n' % (
                          default_tag_file))

    if args.tax:
        if args.tax == '/dev/null':
            sys.stderr.write('[-] Using no taxonomy\n')
        else:
            sys.stderr.write('[-] Using taxonomy in %s\n' % (
                              args.tax))
    else:
        sys.stderr.write('[-] Using default taxonomy in %s\n' % (
                          default_tax_file))

    if args.exp:
        if args.exp == '/dev/null':
            sys.stderr.write('[-] Using no expansion tags\n')
        else:
            sys.stderr.write('[-] Using expansion tags in %s\n' % (
                              args.exp))
    else:
        sys.stderr.write('[-] Using default expansion tags in %s\n' % (
                          default_exp_file))

    main(args)
