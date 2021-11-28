
gawk 'match($0,"^track name=(.{1,50})",m){file=gensub(/[ \"\/\\]/,"_","g",gensub(/\"/,"","g",m[1]))".bed";}{print > file;}' $1

