HTML_FILES := $(patsubst %.md, %.html, $(shell find . -name '*.md' \! -name 'navbar.md'))

all: $(HTML_FILES)

%.html: %.md navbar.html template.html
	pandoc '$<' -o '$@' --smart --template template.html --include-before navbar.html

navbar.html: navbar.md
	pandoc '$<' -o '$@' --smart

.PHONY: clean
clean:
	$(RM) $(HTML_FILES) navbar.html

