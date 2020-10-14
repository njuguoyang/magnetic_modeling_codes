pro check_indexing

    DIR='./'
    OLD='~/DAVE-DAVE4VM_OLD/SHOOTOUT/'
    file='compare.fts'
    work=mrdfits(DIR+file,1)
    base=mrdfits(OLD+file,1)
    stop

end
