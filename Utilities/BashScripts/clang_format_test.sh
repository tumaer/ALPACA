#!/bin/bash

# Check wheter clang-format executable is user defined
terminal_input=$1
clang_format_exec="clang-format"
if [ "$terminal_input" ]
then
  clang_format_exec=$terminal_input
fi

# Applies clang-format

# check that we are in a clean state in order to prevent accidential changes
if [ ! -z "$(git status --untracked-files=no  --porcelain)" ]; then
  echo "Script must be applied on a clean git state"
  exit 1
fi

echo
echo "Checking formatting using the following clang-format version:"
$clang_format_exec --version
echo

# perform clang-format on all cpp-files
find src/ -name '*.h' -or -name '*.hpp' -or -name '*.cpp' | xargs $clang_format_exec -i -style=file
find test/ -name '*.h' -or -name '*.hpp' -or -name '*.cpp' | xargs $clang_format_exec -i -style=file
# check if something was modified
notcorrectlist=`git status --porcelain | grep '^ M' | cut -c4-`
# if nothing changed ok
if [[ -z $notcorrectlist ]]; then
  # send a negative message to gitlab
  echo "Excellent. Very good formatting!"
  exit 0;
else
  echo "The following files have clang-format problems:"
  git diff --stat $notcorrectlist
  echo "Please run"
  echo
  echo "find src/ -name '*.h' -or -name '*.hpp' -or -name '*.cpp' | xargs clang-format -i -style=file"
  echo
  echo "and"
  echo
  echo "find test/ -name '*.h' -or -name '*.hpp' -or -name '*.cpp' | xargs clang-format -i -style=file"
  echo
  echo "to solve the issue."
  # cleanup changes in git
  git reset HEAD --hard
fi

exit 1
