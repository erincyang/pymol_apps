def annotate_good():
    """
    Will save the the objectname(s) in sele to 
    annotate.txt file
    """
    myspace = {'currobj': []}
    cmd.iterate('(sele)', 'currobj.append(model)', space=myspace)
    objs = set(myspace['currobj'])
    with open('annotate.txt','a') as f_out:
        for obj in objs:
            f_out.write(obj + ' 1\n')
            print("annotated good"+obj)

def annotate_bad():
    """
    Will save the the objectname(s) in sele to 
    annotate.txt file
    """
    myspace = {'currobj': []}
    cmd.iterate('(sele)', 'currobj.append(model)', space=myspace)
    objs = set(myspace['currobj'])
    with open('annotate.txt','a') as f_out:
        for obj in objs:
            f_out.write(obj + ' 0\n')
            print("annotated bad"+obj)

cmd.extend("annotate_good", annotate_good)
cmd.extend("annotate_bad", annotate_bad)

