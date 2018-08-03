#!/usr/bin/env bash

echo "#!/usr/bin/env bash" >setKM.sh
if [ -f ~/.kmasker_user.config ]; then
  grep "^export " ~/.kmasker_user.config >>setKM.sh
  chmod 700 setKM.sh
  ./setKM.sh
  rm setKM.sh
fi
