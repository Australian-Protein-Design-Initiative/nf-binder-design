process UNIQUE_ID {
    publishDir path: "${params.outdir}", pattern: "unique_id.txt", mode: 'copy'

    output:
    path 'unique_id.txt', emit: id_file

    script:
    """
    #!/usr/bin/env python3
    import secrets
    
    # Base58 alphabet (excludes ambiguous characters)
    ALPHABET = '123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz'
    
    # Generate 7 cryptographically secure random characters
    unique_id = ''.join(secrets.choice(ALPHABET) for _ in range(7))
    
    with open('unique_id.txt', 'w') as f:
        f.write(unique_id)
    """
}