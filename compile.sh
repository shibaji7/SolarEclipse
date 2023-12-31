isort -rc -sl py/
autoflake --remove-all-unused-imports -i -r py/
isort -rc -m 3 py/
black py

isort -rc -sl pyrt/
autoflake --remove-all-unused-imports -i -r pyrt/
isort -rc -m 3 pyrt/
black pyrt
