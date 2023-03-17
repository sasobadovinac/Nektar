SET(DEVGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(DEVGUIDE ${CMAKE_BINARY_DIR}/docs/developer-guide)

FILE(MAKE_DIRECTORY ${DEVGUIDE}/html)

FIND_PROGRAM(PDFLATEX pdflatex)
MARK_AS_ADVANCED(PDFLATEX)
FIND_PROGRAM(BIBTEX bibtex)
MARK_AS_ADVANCED(BIBTEX)
FIND_PROGRAM(MAKEINDEX makeindex)
MARK_AS_ADVANCED(MAKEINDEX)
FIND_PROGRAM(HTLATEX htlatex)
MARK_AS_ADVANCED(HTLATEX)

ADD_CUSTOM_TARGET(developer-guide-html
    TEXINPUTS=${DEVGUIDESRC}//:
    ${HTLATEX} ${DEVGUIDESRC}/developer-guide.tex 
    "${DEVGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    COMMAND BIBINPUTS=${DEVGUIDESRC} ${BIBTEX} developer-guide.aux
    COMMAND TEXINPUTS=${DEVGUIDESRC}//:
    ${HTLATEX} ${DEVGUIDESRC}/developer-guide.tex 
    "${DEVGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    WORKING_DIRECTORY ${DEVGUIDE}/html
)

# If tex4ht successful, create img dir and copy images across
FILE(GLOB_RECURSE imgfiles RELATIVE ${DEVGUIDESRC}
     "${DEVGUIDESRC}/*.png" "${DEVGUIDESRC}/*.jpg")
FOREACH(img ${imgfiles})
    GET_FILENAME_COMPONENT(DIRNAME ${img} DIRECTORY)
    ADD_CUSTOM_COMMAND(TARGET developer-guide-html POST_BUILD
        COMMENT "Copying ${img} to ${DEVGUIDE}/html/${img}"
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DEVGUIDE}/html/${DIRNAME}
        COMMAND ${CMAKE_COMMAND} -E copy ${img} ${DEVGUIDE}/html/${img}
        WORKING_DIRECTORY ${DEVGUIDESRC})
ENDFOREACH()

FILE(GLOB_RECURSE pdffiles RELATIVE ${DEVGUIDESRC} "*/img/*.pdf")
FIND_PROGRAM(CONVERT convert)
FOREACH(pdf ${pdffiles})
    GET_FILENAME_COMPONENT(BASE ${pdf} NAME_WE)
    GET_FILENAME_COMPONENT(DIR  ${pdf} DIRECTORY)
    # Set custom output location for PDFs in $DEVGUIDE/html/introduction/img/
    IF(${DIR} STREQUAL "introduction/img")
        SET(DIR "img")
    ENDIF()
    ADD_CUSTOM_COMMAND(TARGET developer-guide-html POST_BUILD
        COMMENT "Converting ${pdf} to ${DEVGUIDE}/html/${DIR}/${BASE}.png"
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DEVGUIDE}/html/${DIR}
        COMMAND ${CONVERT} ${pdf} ${DEVGUIDE}/html/${DIR}/${BASE}.png
        WORKING_DIRECTORY ${DEVGUIDESRC})
ENDFOREACH()

ADD_CUSTOM_TARGET(developer-guide-pdf
    ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developer-guide.tex
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${BIBTEX} ${DEVGUIDE}/developer-guide.aux
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${MAKEINDEX} ${DEVGUIDE}/developer-guide.idx
    COMMAND ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developer-guide.tex
    COMMAND ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developer-guide.tex
    WORKING_DIRECTORY ${DEVGUIDESRC}
)

# Install output, if it is generated.
INSTALL(CODE "IF (EXISTS ${DEVGUIDE}/html)
    FILE(GLOB DGHTML ${DEVGUIDE}/html/*)
    LIST(LENGTH DGHTML DGHTML_LEN)
    IF(DGHTML_LEN GREATER 0)
        FILE(INSTALL ${DEVGUIDE}/html
             DESTINATION \${CMAKE_INSTALL_PREFIX}/${NEKTAR_DOC_DIR}/developer-guide)
    ENDIF()
ENDIF()
IF (EXISTS ${DEVGUIDE}/developer-guide.pdf)
    FILE(INSTALL ${DEVGUIDE}/developer-guide.pdf
         DESTINATION \${CMAKE_INSTALL_PREFIX}/${NEKTAR_DOC_DIR}/developer-guide)
ENDIF()")
